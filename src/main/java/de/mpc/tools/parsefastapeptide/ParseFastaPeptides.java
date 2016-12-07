package de.mpc.tools.parsefastapeptide;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;


public class ParseFastaPeptides {

    /** the name of the fasta file */
    private String fastaFileName;

    /** mapping from the peptide to the accessions */
    private Map<String, Set<String>> peptideAccessionMap;

    /** mapping from the peptide to the number of all occurrences (including double occurrences in one protein) */
    private Map<String, Integer> peptideAllOccurrences;

    /** the minimal length of an output peptide */
    private int minLength;

    /** the maximal length of an output peptide */
    private int maxLength;

    /** the number of allowed missed cleavages */
    private int missedCleavages;

    /** used protein digester(s) */
    private ProteinDigester[] digester;


    public ParseFastaPeptides(String fileName) {
        this(fileName, Enzyme.TRYPSIN.toString(), 7, 45, 0);
    }


    public ParseFastaPeptides(String fileName, String enzymeName, int minLength, int maxLength, int missedCleavages) {
        this(fileName, new String[]{enzymeName}, minLength, maxLength, missedCleavages);
    }

    public ParseFastaPeptides(String fileName, String[] enzymeNames, int minLength, int maxLength, int missedCleavages) {

        this.fastaFileName = fileName;
        this.minLength = minLength;
        this.maxLength = maxLength;
        this.missedCleavages = missedCleavages;

        this.digester = new ProteinDigester[enzymeNames.length];
        int c = 0;
        for (String enzymeName : enzymeNames) {
            Enzyme enzyme = Enzyme.valueOf(enzymeName);
            // don't use max length restriction when digesting with multiple enzymes, these will be filtered after digestion
            int max = (enzymeNames.length > 1) ? 0 : this.maxLength;
            this.digester[c++] = new ProteinDigester(enzyme, this.minLength, max, this.missedCleavages);
        }

        if ((this.digester.length > 1) && (missedCleavages > 0)) {
            System.err.println("WARNING: multiple enzymes and missedCleavages > 0 will not work correctly!");
        }
    }


    public void parseFastaFile() throws IOException, DigestException {
        FileInputStream fileStream = new FileInputStream(fastaFileName);

        DataInputStream in = new DataInputStream(fileStream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));

        peptideAccessionMap = new HashMap<String, Set<String>>(10000);
        peptideAllOccurrences = new HashMap<String, Integer>(10000);

        String strLine = null;
        StringBuilder proteinSequence = null;
        String accession = null;

        // first round: digest with the first enzyme
        System.out.println("digesting with " + digester[0].getEnzyme().toString());
        while ((strLine = br.readLine()) != null) {
            if (strLine.startsWith(">")) {
                if ((accession != null) && (proteinSequence != null) && (proteinSequence.length() > 0)) {
                    for (String peptide : digester[0].digest(proteinSequence.toString())) {
                        Set<String> accSet = peptideAccessionMap.get(peptide);
                        if (accSet == null) {
                            accSet = new HashSet<String>();
                            peptideAccessionMap.put(peptide, accSet);
                            peptideAllOccurrences.put(peptide, 0);
                        }
                        accSet.add(accession);
                        peptideAllOccurrences.put(peptide, peptideAllOccurrences.get(peptide) + 1);
                    }
                }

                // start of a new protein
                accession = strLine.split(" ", 2)[0].substring(1);

                proteinSequence = new StringBuilder();
            } else {
                // just reading in the protein sequence
                proteinSequence.append(strLine.trim());
            }
        }
        br.close();
        System.out.println("digestion with " + digester[0].getEnzyme().toString() + " done, " + peptideAccessionMap.size() + " peptides");

        // digest with further enzymes, if any are given
        for (int i=1; i < digester.length; i++) {
            Map<String, Set<String>> roundPeptideAccessionMap = new HashMap<String, Set<String>>(peptideAccessionMap.size());

            System.out.println("digesting with " + digester[i].getEnzyme().toString());
            int pepCount = 0;
            for (Map.Entry<String, Set<String>> peptideIt : peptideAccessionMap.entrySet()) {
                for (String peptide : digester[i].digest(peptideIt.getKey())) {
                    Set<String> accSet = roundPeptideAccessionMap.get(peptide);
                    if (accSet == null) {
                        accSet = new HashSet<String>();
                        roundPeptideAccessionMap.put(peptide, accSet);
                    }
                    accSet.addAll(peptideIt.getValue());

                    if (!peptideAllOccurrences.containsKey(peptide)) {
                        peptideAllOccurrences.put(peptide, 0);
                    }

                    peptideAllOccurrences.put(peptide, peptideAllOccurrences.get(peptide) + peptideIt.getValue().size());
                }

                // reset the counts (even uncut peptides are counted again later)
                peptideAllOccurrences.put(peptideIt.getKey(), peptideAllOccurrences.get(peptideIt.getKey()) - peptideIt.getValue().size());

                pepCount++;
                if (pepCount % 500000 == 0) {
                    System.out.println("\t" + pepCount + " peptides digested with " + digester[i].getEnzyme().toString());
                }
            }

            peptideAccessionMap = roundPeptideAccessionMap;
            System.out.println("digestion with " + digester[i].getEnzyme().toString() + " done, " + peptideAccessionMap.size() + " peptides");
        }

        if ((digester.length > 1) && (maxLength > 0)) {
            System.out.println("removing peptides longer than " + maxLength);
            // filter out the long peptides, if more than one enzyme was used
            Iterator<Map.Entry<String, Set<String>>> mapIt = peptideAccessionMap.entrySet().iterator();
            int remCount = 0;
            while (mapIt.hasNext()) {
                Map.Entry<String, Set<String>> entry = mapIt.next();

                if (entry.getKey().length() > maxLength) {
                    mapIt.remove();
                    remCount++;

                    peptideAllOccurrences.remove(entry.getKey());
                }
            }
            System.out.println("removed " + remCount + " peptides, " + peptideAccessionMap.size() + " remaining");
        }
    }


    public Map<String, Set<String>> getPeptideAccessionMap() {
        return peptideAccessionMap;
    }


    public Integer getPeptideAllOccurences(String peptide) {
        return peptideAllOccurrences.get(peptide);
    }


    public static void main(String[] argv) throws IOException, DigestException {
        ParseFastaPeptides parser = new ParseFastaPeptides(
                "/mnt/data/uniNOBACKUP/FASTAs/18mix_db_plus_contaminants_20081209-edited_nicer_headers.fasta",
                new String[]{Enzyme.TRYPSIN.toString()},
                6,
                0,
                2);

        String outFileName = "/mnt/data/uniNOBACKUP/FASTAs/18mix_db_plus_contaminants_20081209-edited_nicer_headers.fasta.trypsin-2_missed-6_45.txt";
        boolean split = false;
        int splitlength = 1000000;

        parser.parseFastaFile();

        System.out.println("Start writing results to file...");
        FileOutputStream fileStream = new FileOutputStream(outFileName);
        DataOutputStream out = new DataOutputStream(fileStream);
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out));

        String headerline = "peptide\tpeptideLength\t#accessions\t#occurrences\taccessions";
        bw.append(headerline);
        bw.newLine();

        StringBuilder line = null;

        int splitfile = 1;
        int splitLineCount = 0;
        BufferedWriter splitBW = null;
        if (split && (parser.getPeptideAccessionMap().size() > splitlength)) {
            FileOutputStream splitFileStream = new FileOutputStream(outFileName + ".split" + splitfile);
            DataOutputStream splitStream = new DataOutputStream(splitFileStream);
            splitBW = new BufferedWriter(new OutputStreamWriter(splitStream));

            splitBW.append(headerline);
            splitBW.newLine();
        }

        for (Map.Entry<String, Set<String>> mapIt : parser.getPeptideAccessionMap().entrySet()) {
            line = new StringBuilder();

            line.append(mapIt.getKey());
            line.append('\t');
            line.append(mapIt.getKey().length());
            line.append('\t');
            line.append(mapIt.getValue().size());
            line.append('\t');
            line.append(parser.getPeptideAllOccurences(mapIt.getKey()));
            line.append('\t');
            int nrAccs = 0;
            for (String acc : mapIt.getValue()) {
                if (nrAccs++ > 0) {
                    line.append(',');
                }
                line.append(acc);
            }

            bw.append(line);
            bw.newLine();

            if (splitBW != null) {
                splitLineCount++;
                if (splitLineCount % splitlength == 0) {
                    splitBW.close();

                    splitfile++;
                    FileOutputStream splitFileStream = new FileOutputStream(outFileName + ".split" + splitfile);
                    DataOutputStream splitStream = new DataOutputStream(splitFileStream);
                    splitBW = new BufferedWriter(new OutputStreamWriter(splitStream));

                    splitBW.append(headerline);
                    splitBW.newLine();
                }

                splitBW.append(line);
                splitBW.newLine();
            }
        }

        bw.close();

        if (splitBW != null) {
            splitBW.close();
        }

        System.out.println("Results written to file.");
    }
}