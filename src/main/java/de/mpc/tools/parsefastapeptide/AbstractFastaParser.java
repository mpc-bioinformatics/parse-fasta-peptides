package de.mpc.tools.parsefastapeptide;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;


/**
 * Abstract class for FASTA file parsing. Provides the basic parsing functions,
 * and calls {@link #processEntry(String, StringBuilder)} for each sucessfully
 * parsed protein sequence with at least one character.
 *
 * @author julian
 *
 */
public abstract class AbstractFastaParser {

    private String fileName;

    /**
     * Basic constructor setting the filename of the FASTA file
     * @param fileName
     */
    public AbstractFastaParser(String fileName) {
        this.fileName = fileName;
    }


    /**
     * Getter for the file name of the FASTA file
     * @return
     */
    public String getFileName() {
        return fileName;
    }


    /**
     * This function gets the accession and the peptide sequence and processes
     * it further
     */
    public abstract void processEntry(String header, StringBuilder proteinSequence);


    /**
     * Parses a FASTA file and calls processEntry for each protein entry with a
     * sequence longer than 0 characters.
     *
     * @return
     * @throws IOException
     */
    public int parseFastaFile() throws IOException {
        int entryCount = 0;

        LineIterator it = FileUtils.lineIterator(new File(getFileName()), "UTF-8");
        StringBuilder proteinSequence = null;
        String header = null;

        while (it.hasNext()) {
            String strLine = it.nextLine();

            if (strLine.startsWith(">")) {
                checkAndProcess(header, proteinSequence);

                // start of a new entry
                header = strLine.substring(1).trim();
                proteinSequence = new StringBuilder();
                entryCount++;
            } else if (proteinSequence != null) {
                // just reading in the protein sequence
                proteinSequence.append(strLine.trim());
            }
        }

        // process the last entry
        checkAndProcess(header, proteinSequence);

        it.close();
        return entryCount;
    }


    /**
     * Check whether a protein sequence is given and if it is longer than 0.
     * call {@link #processEntry(String, StringBuilder)}, if it is.
     *
     * @param header the FASTA header text, can be null
     * @param proteinSequence the protein sequence to process
     */
    private void checkAndProcess(String header, StringBuilder proteinSequence) {
        if ((proteinSequence != null) && (proteinSequence.length() > 0)) {
            processEntry(header, proteinSequence);
        }
    }
}
