package de.mpc.tools.parsefastapeptide.neo4j;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.graphdb.schema.Schema;
import org.neo4j.unsafe.batchinsert.BatchInserter;
import org.neo4j.unsafe.batchinsert.BatchInserters;

import de.mpc.tools.parsefastapeptide.AbstractFastaParser;
import de.mpc.tools.parsefastapeptide.DigestException;
import de.mpc.tools.parsefastapeptide.Enzyme;
import de.mpc.tools.parsefastapeptide.ProteinDigester;
import uk.ac.ebi.pride.utilities.mol.MoleculeUtilities;
import uk.ac.ebi.pride.utilities.mol.NeutralLoss;
import uk.ac.ebi.pride.utilities.mol.NuclearParticle;


public class ParseToNeo4J extends AbstractFastaParser {

    private static final Logger LOGGER = LogManager.getLogger("ParseToNeo4J");


    /** the used enzyme */
    private ProteinDigester enzyme;

    /** path to the Neo4j database */
    private String dbPath;

    /** minimal peptide length */
    private int minLength = 6;

    /** maximal peptide length */
    private int maxLength = 45;

    /** number of missed cleavages allowed by the enzyme */
    private int missedCleavages = 2;

    /** charges for the ions */
    private Integer[] charges = new Integer[]{2, 3};

    /** considered fixed modifications */
    private Map<Character, Double> fixedModifications;

    /** considered variable modifications */
    private Map<Character, Double> variableModifications;


    /** whether the unmodified ion should be stored as well, if fixed modifications are given */
    private boolean encodeUnmodified = true;

    /** counter for the processed accessions */
    private long processedAccessions;



    /** mapping from peptide sequences to the position of the accessions in this batch-round only */
    private Map<String, Set<Integer>> peptidesInBatch;

    /** the accessions in this batch */
    private List<String> accessionsInBatch;

    /** counter for added peptides*/
    private long addedPeptides;

    /** maximal number of peptides before issueing an insertion */
    private long maxPeptidesBeforeInsert;

    private boolean pepSequenceIndexCreated;


    // constants
    private static final Label LABEL_ACCESSION = Label.label("accession");
    private static final Label LABEL_PEPTIDE = Label.label("peptide");
    private static final Label LABEL_ION = Label.label("ion");

    private static final String PROPERTY_SEQUENCE = "sequence";
    private static final String PROPERTY_ACCESSION = "accession";
    private static final String PROPERTY_DESCRIPTION = "description";
    private static final String PROPERTY_LENGTH = "length";
    private static final String PROPERTY_CHARGE = "charge";
    private static final String PROPERTY_MODIFICATION_FIXED = "modification_fixed";
    private static final String PROPERTY_MODIFICATION_VARIABLE = "modification_variable";
    private static final String PROPERTY_MASS_TO_CHARGE = "mass_to_charge";

    private static final CharSequence FIXED_MODIFICATION_SPLITTER = "---";


    // user: neo4j
    // passwd: graph


    /**
     * Constructor
     *
     * @param fileName
     */
    public ParseToNeo4J(String fileName, String dbPath) {
        super(fileName);

        this.dbPath = dbPath;

        // TODO: make these settings adjustable
        // >>>>>>> settings from here
        enzyme = new ProteinDigester(Enzyme.TRYPSIN, minLength, maxLength, missedCleavages);
        //enzyme = new ProteinDigester(Enzyme.CUTALL, minLength, maxLength, 45);

        fixedModifications = new HashMap<>();
        fixedModifications.put('C', 57.021464);

        variableModifications = new HashMap<>();
        variableModifications.put('M', 15.994915);

        maxPeptidesBeforeInsert = 5000000;
        // <<<<<<< settings up to here

        processedAccessions = 0;
        addedPeptides = 0;
        pepSequenceIndexCreated = false;

        // caching of peptides for the insertion round
        peptidesInBatch = new HashMap<>();
        accessionsInBatch = new ArrayList<>();

        initializeGraphDB();
    }


    /**
     * sets the maximal number of parsed peptides before inserting into the DB.
     *
     * @param maxPeps
     */
    public void setMaxPeptidesBeforeInsert(long maxPeps) {
        this.maxPeptidesBeforeInsert = maxPeps;
        LOGGER.info("set number of peptides before inserting into DB to {}", maxPeptidesBeforeInsert);
    }


    /**
     * Shuts down the database connection
     */
    public void shutdown() {
        LOGGER.info("Cleaning up temp files.");

        // nothing on file level to do anymore

        LOGGER.info("Cleanup done.");
    }


    /**
     * Registers a shutdown hook for the Neo4j instance so that it shuts down
     * nicely when the VM exits (even if you "Ctrl-C" the running application).
     *
     * @param graphDb
     */
    private static void registerShutdownHook( final GraphDatabaseService graphDb ) {
        Runtime.getRuntime().addShutdownHook( new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        } );
    }


    /**
     * Connect to the database
     */
    private void initializeGraphDB() {
        // create store for the graph DB
        File dbPathFile = new File(dbPath);
        if (!dbPathFile.exists()) {
            try {
                Files.createDirectory(dbPathFile.toPath());
            } catch (IOException ex) {
                LOGGER.error("Cannot create directory '{}'", dbPath, ex);
            }
        }
    }


    @Override
    public int parseFastaFile() throws IOException {
        int parsedEntries = super.parseFastaFile();

        // insert the last peptides
        insertPeptidesOfBatchIntoDB();
        LOGGER.info("Added {} peptides to the DB.", addedPeptides);

        // add the indizes
        createIndizes();

        return parsedEntries;
    }


    @Override
    public void processEntry(String header, StringBuilder proteinSequence) {
        Integer accessionID = accessionsInBatch.size();
        accessionsInBatch.add(header);

        try {
            // digest the sequence and cache the peptides
            enzyme.digest(proteinSequence.toString())
                    .forEach(peptide -> {
                        if (MoleculeUtilities.isAminoAcidSequence(peptide)) {
                            // save accession Id for the peptide
                            peptidesInBatch.computeIfAbsent(peptide, pep -> new HashSet<Integer>()).add(accessionID);
                        } else {
                            LOGGER.error("Could not add peptide for '{}', this is considered to be no peptide sequence: '{}'", header, peptide);
                        }
                    });
        } catch (DigestException e) {
            LOGGER.error("error digesting sequence: {}", proteinSequence, e);
        }

        processedAccessions++;
        if (processedAccessions % 10000 == 0) {
            LOGGER.info("processed {} entries ({} peptides in batch)...", processedAccessions, peptidesInBatch.size());
        }

        if (peptidesInBatch.size() >= maxPeptidesBeforeInsert) {
            insertPeptidesOfBatchIntoDB();
        }
    }


    /**
     * Processes the peptides of the batch into the graphDB
     */
    private void insertPeptidesOfBatchIntoDB() {
        LOGGER.info("insert batch called with {} accessions and {} peptides", processedAccessions, peptidesInBatch.size());

        Map<String, Long> pepMapSequenceToIDs = getPeptideIDsInGraph(peptidesInBatch.keySet());

        try {
            BatchInserter batchInserter = BatchInserters.inserter( new File(dbPath) );

            // add accessions to graph and map to the accession IDs in the graph
            List<Long> accessionIds = accessionsInBatch.stream()
                    .map(header -> addAccessionToGraph(header, batchInserter))
                    .collect(Collectors.toList());

            int count = 0;
            for (Map.Entry<String, Set<Integer>> pepInBatch : peptidesInBatch.entrySet()) {
                // get peptide's ID or create it
                Long pepNodeId = pepMapSequenceToIDs.get(pepInBatch.getKey());
                if (pepNodeId == null) {
                    pepNodeId = insertPeptideInDB(pepInBatch.getKey(), batchInserter);
                }

                // connect to accessions
                for (Integer acc : pepInBatch.getValue()) {
                    Long accNodeId = accessionIds.get(acc);
                    batchInserter.createRelationship(accNodeId, pepNodeId, DigestedRelTypes.BELONGS_TO, null);
                }

                if (++count % (maxPeptidesBeforeInsert / 10) == 0) {
                    LOGGER.info("added {} of {} peptides", count, peptidesInBatch.size());
                }
            }

            LOGGER.info("Closing batch inserter.");
            batchInserter.shutdown();
            LOGGER.info("Batch inserter closed.");
        } catch (IOException e) {
            LOGGER.error("error while inserting batch of peptides", e);
        }

        // clearing the batch entries
        peptidesInBatch.clear();
        accessionsInBatch.clear();
    }


    /**
     * Gets the peptide Ids for the given peptides, which are already in teh graph.
     * @param peptides
     * @return mapping of the given peptides, which are already in the graph
     */
    private Map<String, Long> getPeptideIDsInGraph(Collection<String> peptides) {
        LOGGER.info("Getting IDs of batch's peptides, that are already in the DB");

        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase( new File(dbPath) );
        registerShutdownHook(graphDb);

        // create an index on the peptides
        if (!pepSequenceIndexCreated) {
            createPeptideSequenceIndex(graphDb);
        }

        Map<String, Long> pepIdMap = new HashMap<>(peptides.size());

        try (Transaction tx = graphDb.beginTx()) {
            for (String peptide : peptides) {
                Node pepNode = graphDb.findNode(LABEL_PEPTIDE, PROPERTY_SEQUENCE, peptide);
                if (pepNode != null) {
                    pepIdMap.put(peptide, pepNode.getId());
                }
            }

            tx.success();
        }

        LOGGER.info("found {} peptides already in the DB", pepIdMap.size());

        //dropPeptideIndizes(graphDb);

        graphDb.shutdown();
        return pepIdMap;
    }


    /**
     * Adds the given accession into the graph using the {@link BatchInserter}
     *
     * @param header
     * @param batchInserter
     * @return
     */
    private long addAccessionToGraph(String header, BatchInserter batchInserter) {
        String[] splitHeader = header.split("\\s", 2);
        Map<String, Object> accProperties = new HashMap<>();
        accProperties.put(PROPERTY_ACCESSION, splitHeader[0]);
        if (splitHeader.length > 1) {
            accProperties.put(PROPERTY_DESCRIPTION, splitHeader[1]);
        } else {
            accProperties.put(PROPERTY_DESCRIPTION, splitHeader[0]);
        }

        //accProperties.put(PROPERTY_SEQUENCE, proteinSequence );
        return batchInserter.createNode(accProperties, LABEL_ACCESSION);
    }


    /**
     * Inserts the peptide with the given sequence in the DB, together with its ions.
     *
     * @param peptide th epeptide sequence
     * @return the peptide node id
     */
    private long insertPeptideInDB(String peptide, BatchInserter batchInserter) {
        Map<String, Object> pepProperties = new HashMap<>(2);
        pepProperties.put(PROPERTY_SEQUENCE, peptide);
        pepProperties.put(PROPERTY_LENGTH, peptide.length());
        long pepID = batchInserter.createNode(pepProperties, LABEL_PEPTIDE);

        addIonsForPeptide(pepID, peptide, batchInserter);

        addedPeptides++;
        return pepID;
    }


    /**
     * Adds the ions created by the given peptide to the database. All possible charges and modifications are iterated.
     *
     * @param pepNodeId
     * @param sequence
     */
    private void addIonsForPeptide(long pepNodeId, String sequence, BatchInserter batchInserter) {
        Map<String, Double> modificationMasses = getPossibleModificationMasses(sequence);

        Double enzymeAddedMass = NeutralLoss.WATER_LOSS.getMonoMass();
        Double unmodifiedMass = MoleculeUtilities.calculateTheoreticalMass(sequence, enzymeAddedMass);

        Map<String, Object> ionProperties = new HashMap<>(1);
        Map<String, Object> ionisationProperties = new HashMap<>(3);

        modificationMasses.forEach((encMods, massShift) ->  {
            double theoreticalMass = unmodifiedMass + massShift;

            for (Integer charge : charges) {
                double theoreticalMz = (theoreticalMass + charge * NuclearParticle.PROTON.getMonoMass()) / charge;

                ionProperties.clear();
                ionProperties.put(PROPERTY_MASS_TO_CHARGE, theoreticalMz);
                long ionID = batchInserter.createNode(ionProperties, LABEL_ION);

                ionisationProperties.clear();
                ionisationProperties.put(PROPERTY_CHARGE, charge);
                String[] splittedModifications = splitEncodedModifications(encMods, massShift);
                ionisationProperties.put(PROPERTY_MODIFICATION_FIXED, splittedModifications[0]);
                ionisationProperties.put(PROPERTY_MODIFICATION_VARIABLE, splittedModifications[1]);

                batchInserter.createRelationship(pepNodeId, ionID, DigestedRelTypes.BELONGS_TO, ionisationProperties);
            }
        });
    }


    /**
     * Split the encoded modifications into the fixed and variable part.
     *
     * @param encodedModifications
     * @param totalMassShift
     * @return
     */
    private String[] splitEncodedModifications(String encodedModifications, double totalMassShift) {
        String fixedModifcations = "";
        String variableModifcations = "";

        if ((totalMassShift > 0) || (totalMassShift < 0)) {
            if (encodedModifications.contains(FIXED_MODIFICATION_SPLITTER)) {
                String[] splittedModifcations = encodedModifications.split(FIXED_MODIFICATION_SPLITTER.toString());
                fixedModifcations = splittedModifcations[0];
                if (splittedModifcations.length > 1) {
                    variableModifcations = splittedModifcations[1];
                }
            } else  {
                fixedModifcations = "";
                variableModifcations = encodedModifications;
            }
        }

        return new String[]{fixedModifcations, variableModifcations};
    }


    /**
     * Creates a mapping for all possible modifications of the given sequence.
     *
     * @param sequence
     * @return mapping from encoded modifications to the mass shift
     */
    private Map<String, Double> getPossibleModificationMasses(String sequence) {
        Map<String, Double> modificationMasses = new HashMap<>();
        modificationMasses.put("", 0.0);

        double fixedMass = 0;
        Set<Character> foundFixed = new HashSet<>(fixedModifications.size());
        if (!fixedModifications.isEmpty() || !variableModifications.isEmpty()) {
            char[] seqChar = sequence.toCharArray();
            for (int idx=0; idx < seqChar.length; idx++) {
                if (fixedModifications.keySet().contains(seqChar[idx])) {
                    fixedMass += fixedModifications.get(seqChar[idx]);
                    foundFixed.add(seqChar[idx]);
                }

                if (variableModifications.keySet().contains(seqChar[idx])) {
                    double mass = variableModifications.get(seqChar[idx]);

                    StringBuilder variableString = new StringBuilder();
                    variableString.append(seqChar[idx]);
                    variableString.append(idx+1);
                    variableString.append('[');
                    variableString.append(mass);
                    variableString.append(']');

                    modificationMasses = addVariableModificationToMap(mass,
                            variableString.toString(), modificationMasses);
                }
            }
        }

        // now add the fixed modifications
        if (!foundFixed.isEmpty()) {
            StringBuilder fixedSb = new StringBuilder();
            for (Character fixed : foundFixed) {
                fixedSb.append(fixed);
                fixedSb.append('[');
                fixedSb.append(fixedModifications.get(fixed));
                fixedSb.append(']');
            }
            String fixedString = fixedSb.toString();

            Map<String, Double> iterMasses = modificationMasses;
            modificationMasses = new HashMap<>(iterMasses.size());
            for (Map.Entry<String, Double> modIt : iterMasses.entrySet()) {
                modificationMasses.put(fixedString + FIXED_MODIFICATION_SPLITTER + modIt.getKey(),
                        fixedMass + modIt.getValue());

                if (encodeUnmodified) {
                    modificationMasses.put(modIt.getKey(), modIt.getValue());
                }
            }
        }

        return modificationMasses;
    }


    /**
     * Adds to the given map of encoded modifications a new variable
     * modification with the given mass and encoding. Effectively, the map size
     * will be doubled by this.
     *
     * @param mass
     * @param modEncode
     * @param modificationMasses
     * @return
     */
    private Map<String, Double> addVariableModificationToMap(double mass, String modEncode,
            Map<String, Double> modificationMasses) {
        Map<String, Double> varModMasses = new HashMap<>(modificationMasses.size()*2);

        for (Map.Entry<String, Double> modIt : modificationMasses.entrySet()) {
            varModMasses.put(modIt.getKey() + modEncode, modIt.getValue() + mass);
        }

        varModMasses.putAll(modificationMasses);
        return varModMasses;
    }


    /**
     * Creates the index on the peptide sequence
     */
    private void createPeptideSequenceIndex(GraphDatabaseService graphDb) {
        LOGGER.info("Creating peptide sequence index");

        try (Transaction tx = graphDb.beginTx()) {
            Schema schema = graphDb.schema();

            schema.indexFor(LABEL_PEPTIDE)
                    .on(PROPERTY_SEQUENCE)
                    .create();

            tx.success();
        }

        try (Transaction tx = graphDb.beginTx()) {
            Schema schema = graphDb.schema();
            schema.awaitIndexesOnline(30, TimeUnit.DAYS);
            tx.success();
        }

        pepSequenceIndexCreated = true;

        LOGGER.info("Peptide Sequence index created.");
    }


    /**
     * Creates the indizes for the graphDB
     */
    private void createIndizes() {
        LOGGER.info("Creating indizes.");

        GraphDatabaseService graphDb = new GraphDatabaseFactory().newEmbeddedDatabase( new File(dbPath) );
        registerShutdownHook(graphDb);

        try (Transaction tx = graphDb.beginTx()) {
            Schema schema = graphDb.schema();

            if (!pepSequenceIndexCreated) {
                schema.indexFor(LABEL_PEPTIDE)
                        .on(PROPERTY_SEQUENCE)
                        .create();
            }

            schema.indexFor(LABEL_ACCESSION)
                    .on(PROPERTY_ACCESSION)
                    .create();

            schema.indexFor(LABEL_ION)
                    .on(PROPERTY_MASS_TO_CHARGE)
                    .create();

            tx.success();
        }

        try (Transaction tx = graphDb.beginTx()) {
            Schema schema = graphDb.schema();
            schema.awaitIndexesOnline(30, TimeUnit.DAYS);
            tx.success();
        }
        LOGGER.info("Indizes created.");

        LOGGER.info("Shutting down graphDB connection after indizes' creation.");
        graphDb.shutdown();
        LOGGER.info("Shutted down.");
    }


    public static void main(String[] argv) {
        //String filePath = "/mnt/data/uniNOBACKUP/FASTAs/cRAP-contaminants-better_parseable_headers-20120229.fasta";
        //String filePath = "/mnt/data/uniNOBACKUP/FASTAs/EcoProt_lib_20160518.fasta";
        //String filePath = "/mnt/data/uniNOBACKUP/FASTAs/uniprot_complete-homo_sapiens_20170627.fasta";
        //String filePath = "/mnt/data/uniNOBACKUP/FASTAs/shared/2015_11/uniprot_sprot.fasta";

        String filePath;
        String dbPath = null;

        if (argv.length > 0) {
            filePath = argv[0];
        } else {
            filePath = "/dev/null";
        }

        if (argv.length > 1) {
            dbPath = argv[1];
        } else {
            try {
                Path tmpPath = Files.createTempDirectory("parserGraphDB-");
                dbPath = tmpPath.toAbsolutePath().toString();
            } catch (IOException ex) {
                LOGGER.error("Could not create a tmp folder for the DB", ex);
            }
        }

        ParseToNeo4J parser = new ParseToNeo4J(filePath, dbPath);

        if (argv.length > 2) {
            long maxPeps;
            try {
                maxPeps = Long.parseLong(argv[2]);
            } catch (NumberFormatException e) {
                maxPeps = 5000000;
            }
            parser.setMaxPeptidesBeforeInsert(maxPeps);
        }

        LOGGER.info("start parsing");
        try {
            parser.parseFastaFile();
            LOGGER.info("parsing done");
        } catch (IOException e) {
            LOGGER.error(e);
        }

        parser.shutdown();
    }

    /*

neo4j snippets:

shared peptides:
MATCH (a:accession)-[r:BELONGS_TO]->(p:peptide)
WITH p, count(r) as rel_cnt
WHERE rel_cnt > 1
RETURN count(p)


peptides in a mas ratio and charge
MATCH (p:peptide)-[r:BELONGS_TO]->(i:ion)
WHERE
        i.mass_to_charge >= 299.7
        and i.mass_to_charge <= 300.3
        and r.charge = 3
RETURN count(p)

     */

}
