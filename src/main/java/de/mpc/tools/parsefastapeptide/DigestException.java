package de.mpc.tools.parsefastapeptide;

/**
 * Not very elaborated exception which may occur during the digestion.
 * 
 * @author julian
 *
 */
public class DigestException extends Exception {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	
	public DigestException(String exc) {
		super(exc);
	}
}
