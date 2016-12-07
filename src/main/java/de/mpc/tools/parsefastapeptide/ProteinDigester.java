package de.mpc.tools.parsefastapeptide;

import java.util.List;


/**
 * This will be the protein digester base class (i.e. modeller)
 * 
 * @author julian
 *
 */
public class ProteinDigester {
	/** the used enzyme */
	private Enzyme enzyme;
	
	/** the minimal length of an output peptide */
	private int minLength;
	
	/** the maximal length of an output peptide */
	private int maxLength;
	
	/** the number of allowed missed cleavages */
	private int missedCleavages;
	
	
	public ProteinDigester() {
		this(null, 7, 45, 0);
	}
	
	
	public ProteinDigester(Enzyme enzyme, int minLength, int maxLength, int missedCleavages) {
		this.enzyme = enzyme;
		this.minLength = minLength;
		this.maxLength = maxLength;
		this.missedCleavages = missedCleavages;
	}
	
	
	/**
	 * Setter for the enzyme
	 * 
	 * @param enzyme
	 */
	public void setEnzyme(Enzyme enzyme) {
		this.enzyme = enzyme;
	}
	
	
	/**
	 * Getter for the enzyme
	 * @return
	 */
	public Enzyme getEnzyme() {
		return enzyme;
	}
	
	
	/**
	 * Digest the given protein
	 */
	public List<String> digest(String proteinSequence) throws DigestException {
		// here we finally have some logic 
		if (proteinSequence == null) {
			throw new DigestException("No protein sequence given for digestion.");
		}
		
		if (enzyme == null) {
			throw new DigestException("No enzyme given for digestion.");
		}
		
		return enzyme.digestProtein(proteinSequence.replaceAll("\\s", "").toUpperCase(),
				minLength, maxLength, missedCleavages);
	}
	
	
	/**
	 * Digest the given protein with the given enzyme and returns all peptides
	 * with length between minLength and maxLength
	 * 
	 * @param protein
	 * @param enzyme
	 * @param minLength
	 * @param maxLength
	 * @param missedCleavages
	 * @return either the list of peptides or null if error while digesting
	 */
	public static List<String> digestProtein(String protein, Enzyme enzyme, int minLength, int maxLength, int missedCleavages)
			throws DigestException{
		ProteinDigester digester = new ProteinDigester(enzyme, minLength, maxLength, missedCleavages);
		return digester.digest(protein);
	}
}