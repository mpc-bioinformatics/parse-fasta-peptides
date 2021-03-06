package de.mpc.tools.parsefastapeptide;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class EnzymeTest {
	
	@Test
	public void testParseProtein() {
		
		String protein = "MTEYKLVVVGAAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG"
				+ "REEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNNCDL"
				+ "PSRTVDTKQAQDLARSYGIPFIETSTKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPGC"
				+ "VKIKKCIIM";
		
		int minLength = 6;
		int maxLength = 0;
		
		
		System.out.println(Enzyme.TRYPSIN.digestProtein(protein, minLength, maxLength));
		
		System.out.println("\n\n");
		System.out.println(Enzyme.TRYPSIN.digestProtein(protein, minLength, maxLength, 2));
		
		//assertEquals("trypsin (MS:1001176) should be found in the obo", true, foundTrypsin);
		
		/*
		MTEYK								MTEYK
		LVVVGAAGVGK							LVVVGAAGVGK
											MTEYKLVVVGAAGVGK
		SALTIQLIQNHFVDEYDPTIEDSYR			SALTIQLIQNHFVDEYDPTIEDSYR
											LVVVGAAGVGKSALTIQLIQNHFVDEYDPTIEDSYR
		K
		QVVIDGETCLLDILDTAGR
		EEYSAMR
		DQYMR
		TGEGFLCVFAINNTK
		SFEDIHHYR
		EQIK
		R
		VK
		DSEDVPMVLVGNNCDLPSR
		TVDTK
		QAQDLAR
		SYGIPFIETSTK
		TR
		QR
		VEDAFYTLVR
		EIR
		QYR
		LK
		K
		ISK
		EEK
		TPGCVK
		IK
		K
		CIIM									, , K, SALTIQLIQNHFVDEYDPTIEDSYRK, QVVIDGETCLLDILDTAGR, KQVVIDGETCLLDILDTAGR, EEYSAMR, QVVIDGETCLLDILDTAGREEYSAMR, DQYMR, EEYSAMRDQYMR, TGEGFLCVFAINNTK, DQYMRTGEGFLCVFAINNTK, SFEDIHHYR, TGEGFLCVFAINNTKSFEDIHHYR, EQIK, SFEDIHHYREQIK, R, EQIKR, VK, RVK, DSEDVPMVLVGNNCDLPSR, VKDSEDVPMVLVGNNCDLPSR, TVDTK, DSEDVPMVLVGNNCDLPSRTVDTK, QAQDLAR, TVDTKQAQDLAR, SYGIPFIETSTK, QAQDLARSYGIPFIETSTK, TR, SYGIPFIETSTKTR, QR, TRQR, VEDAFYTLVR, QRVEDAFYTLVR, EIR, VEDAFYTLVREIR, QYR, EIRQYR, LK, QYRLK, K, LKK, ISK, KISK, EEK, ISKEEK, TPGCVK, EEKTPGCVK, IK, TPGCVKIK, K, IKK, CIIM, KCIIM
		*/
	}
}
