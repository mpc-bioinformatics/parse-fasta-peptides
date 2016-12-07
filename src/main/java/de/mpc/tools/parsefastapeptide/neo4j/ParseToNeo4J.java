package de.mpc.tools.parsefastapeptide.neo4j;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.graphdb.schema.IndexDefinition;
import org.neo4j.graphdb.schema.Schema;

import de.mpc.tools.parsefastapeptide.AbstractFastaParser;
import de.mpc.tools.parsefastapeptide.DigestException;
import de.mpc.tools.parsefastapeptide.Enzyme;
import de.mpc.tools.parsefastapeptide.ProteinDigester;
import uk.ac.ebi.pride.utilities.mol.AminoAcid;
import uk.ac.ebi.pride.utilities.mol.AminoAcidSequence;
import uk.ac.ebi.pride.utilities.mol.MoleculeUtilities;
import uk.ac.ebi.pride.utilities.mol.NuclearParticle;


public class ParseToNeo4J extends AbstractFastaParser {

    private static final Logger logger = LogManager.getLogger("ParseToNeo4J");

    private static final String DB_PATH = "/home/julian/opt/neo4j-community-3.0.7/data/databases/graph.db/";


    /** connection to the Neo4J database */
    private GraphDatabaseService graphDb;

    /** the used enzyme */
    private ProteinDigester enzyme;

    /** minimal peptide length */
    private int minLength = 5;

    /** maximal peptide length */
    private int maxLength = 45;

    /** number of missed cleavages allowed by the enzyme */
    private int missedCleavages = 45;

    /** charges for the ions */
    private Integer[] charges = new Integer[]{2, 3};

    /** considered fixed modifications */
    private Map<Character, Double> fixedModifications;

    /** considered variable modifications */
    private Map<Character, Double> variableModifications;

    /** maps from peptide string to peptide */
    private Map<String, Node> pepMap;

    /** whether the unmodified ion should be stored as well, if fixed modifications are given */
    private boolean encodeUnmodified = true;

    /** counter for the processed entries */
    private int processedEntries;

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

    // TODO: include
    // - fixed modifications
    // - variable modifications
    // - charges

    //      accession (sequence, accession, description)
    //          |
    //       peptide  (sequence, length)
    //          |
    //         ion    (charge, fixed mod, variable mod, m/z)


    /**
     * Constructor
     *
     * @param fileName
     */
    public ParseToNeo4J(String fileName) {
        super(fileName);

        // TODO: make these settings adjustable
        enzyme = new ProteinDigester(Enzyme.TRYPSIN, minLength, maxLength, missedCleavages);

        fixedModifications = new HashMap<>();
        fixedModifications.put('C', 57.0);

        variableModifications = new HashMap<>();
        variableModifications.put('M', 18.0);

        pepMap = new HashMap<>(100000);

        processedEntries = 0;

        initializeGraphDB();
    }


    /**
     * Shuts down the database connection
     */
    public void shutdown() {
        if (graphDb != null) {
            graphDb.shutdown();
        }
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
     * Connect to the database and create indizes
     */
    private void initializeGraphDB() {
        // connect to the database
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabase( new File(DB_PATH) );
        registerShutdownHook(graphDb);


        // TODO remove this later, only for convinience
        logger.debug("deleting contents of DB");
        try (Transaction tx = graphDb.beginTx()) {
            for (Relationship rel : graphDb.getAllRelationships()) {
                rel.delete();
            }
            for (Node node : graphDb.getAllNodes()) {
                node.delete();
            }
            tx.success();
        }
        try (Transaction tx = graphDb.beginTx()) {
            for (IndexDefinition idx :  graphDb.schema().getIndexes()) {
                idx.drop();
            }
            tx.success();
        }
        logger.debug("DB cleaned");

        // create indizes
        try (Transaction tx = graphDb.beginTx()) {
            Schema schema = graphDb.schema();

            // indizes on accession
            schema.indexFor(LABEL_ACCESSION)
                    .on(PROPERTY_ACCESSION)
                    .create();
            // indizes on peptide
            schema.indexFor(LABEL_PEPTIDE)
                    .on(PROPERTY_SEQUENCE)
                    .create();
            schema.indexFor(LABEL_PEPTIDE)
                    .on(PROPERTY_LENGTH)
                    .create();
            // indizes on ion
            schema.indexFor(LABEL_ION)
                    .on(PROPERTY_CHARGE)
                    .create();
            schema.indexFor(LABEL_ION)
                    .on(PROPERTY_MODIFICATION_FIXED)
                    .create();
            schema.indexFor(LABEL_ION)
                    .on(PROPERTY_MODIFICATION_VARIABLE)
                    .create();
            schema.indexFor(LABEL_ION)
                    .on(PROPERTY_MASS_TO_CHARGE)
                    .create();

            tx.success();
        }
    }


    @Override
    public void processEntry(String header, StringBuilder proteinSequence) {
        try {
            String[] splitHeader = header.split("\\s", 2);
            processDigestedEntry(splitHeader[0], splitHeader[1],
                    enzyme.digest(proteinSequence.toString()));

            processedEntries++;

            if (processedEntries % 10 == 0) {
                logger.info("processed {} entries...", processedEntries);
            }
        } catch (DigestException e) {
            logger.error("Error while processing '" + header + "'", e);
        }
    }


    /**
     * performs all that is necessary to put the digested entry of the FASTA
     * file into the database.
     *
     * @param accession
     * @param description
     * @param peptides
     */
    private void processDigestedEntry(String accession, String description, List<String> peptides) {
        logger.debug("acc '{}', desc: {}, {}", accession, description, peptides.size());

        try (Transaction tx = graphDb.beginTx()) {
            // create the accession node
            Node accessionNode = graphDb.createNode(LABEL_ACCESSION);
            accessionNode.setProperty(PROPERTY_ACCESSION, accession);
            accessionNode.setProperty(PROPERTY_DESCRIPTION, description);
            //accessionNode.setProperty(PROPERTY_SEQUENCE, sequence);

            for (String pep : peptides) {
                Node pepNode = pepMap.get(pep);

                if (pepNode == null) {
                    pepNode = graphDb.createNode(LABEL_PEPTIDE);
                    pepNode.setProperty(PROPERTY_SEQUENCE, pep);
                    pepNode.setProperty(PROPERTY_LENGTH, pep.length());
                    pepMap.put(pep, pepNode);

                    addIonsForPeptide(pepNode, pep);
                }

                Relationship belongsRel= pepNode.createRelationshipTo(accessionNode,
                        DigestedRelTypes.BELONGS_TO);
                // TODO: maybe add a position relationship
            }

            tx.success();
        }

    }


    /**
     * Adds the ions created by the given paptide to the database. All possible
     * charges and modifications are iterated.
     *
     * @param pepNode
     * @param sequence
     */
    private void addIonsForPeptide(Node pepNode, String sequence) {
        Map<String, Double> modificationMasses = getPossibleModificationMasses(sequence);

        Double unmodifiedMass = MoleculeUtilities.calculateTheoreticalMass(sequence, 0.0);

        logger.debug("ions for {}", sequence);
        for (Map.Entry<String, Double> modIt : modificationMasses.entrySet()) {
            double theoreticalMass = unmodifiedMass + modIt.getValue();

            for (Integer charge : charges) {
                double theoreticalMz = theoreticalMass / charge
                        + NuclearParticle.PROTON.getMonoMass() * charge;

                Node ionNode = graphDb.createNode(LABEL_ION);

                ionNode.setProperty(PROPERTY_CHARGE, charge);
                ionNode.setProperty(PROPERTY_MASS_TO_CHARGE, theoreticalMz);

                String[] splittedModifications = splitEncodedModifications(modIt.getKey(), modIt.getValue());
                ionNode.setProperty(PROPERTY_MODIFICATION_FIXED, splittedModifications[0]);
                ionNode.setProperty(PROPERTY_MODIFICATION_VARIABLE, splittedModifications[1]);

                ionNode.createRelationshipTo(pepNode, DigestedRelTypes.BELONGS_TO);

                logger.debug("\t charge {}, mods: {}, mass {}", charge, modIt, theoreticalMz);
            }
        }
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

        if (totalMassShift != 0) {
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
     * @return mapping from encoded modiifcations to the mass shift
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



    public static void main(String[] argv) {
        //String filePath = "/mnt/data/uniNOBACKUP/FASTAs/cRAP-contaminants-better_parseable_headers-20120229.fasta";
        String filePath = "/mnt/data/uniNOBACKUP/FASTAs/EcoProt_lib_20160518.fasta";

        ParseToNeo4J parser = new ParseToNeo4J(filePath);
        logger.info("start parsing");
        try {
            parser.parseFastaFile();
            logger.info("parsing done");
        } catch (IOException e) {
            logger.error(e);
        }

        parser.shutdown();
    }
}
