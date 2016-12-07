package de.mpc.tools.parsefastapeptide;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;
import java.util.regex.Pattern;

/**
 * This enum defines enzymes for the digestion
 *
 * @author julian
 *
 */
public enum Enzyme {

    /**
     * cut after each amino acid
     */
    CUTALL {
        private String restrictionRules = "(?<=.)";
        private Pattern pattern = Pattern.compile(restrictionRules);

        @Override
        public String getRestrictionRules() {
            return restrictionRules;
        }

        @Override
        protected Pattern getRestrictionPattern() {
            return pattern;
        }
    },

    /**
     * Trypsin cuts like (?<=[KR])(?!P)
     */
    TRYPSIN {
        private String restrictionRules = "(?<=[KR])(?!P)";
        private Pattern pattern = Pattern.compile(restrictionRules);

        @Override
        public String getRestrictionRules() {
            return restrictionRules;
        }

        @Override
        protected Pattern getRestrictionPattern() {
            return pattern;
        }
    },


    /**
     * Chymotrypsin cuts like (?<=[FYWL])(?!P)
     */
    CHYMOTRYPSIN {
        private String restrictionRules = "(?<=[FYWL])(?!P)";
        private Pattern pattern = Pattern.compile(restrictionRules);

        @Override
        public String getRestrictionRules() {
            return restrictionRules;
        }

        @Override
        protected Pattern getRestrictionPattern() {
            return pattern;
        }
    },

    CNBR {
        private String restrictionRules = "(?<=M)";
        private Pattern pattern = Pattern.compile(restrictionRules);

        @Override
        public String getRestrictionRules() {
            return restrictionRules;
        }

        @Override
        protected Pattern getRestrictionPattern() {
            return pattern;
        }
    },

    PROTEINASEK {
        private String restrictionRules = "(?<=[FYWLIAV])";
        private Pattern pattern = Pattern.compile(restrictionRules);

        @Override
        public String getRestrictionRules() {
            return restrictionRules;
        }

        @Override
        protected Pattern getRestrictionPattern() {
            return pattern;
        }
    },

    ;


    /**
     * Returns the pattern of the enzyme's digestion as string
     *
     * @return
     */
    public abstract String getRestrictionRules();


    /**
     * Returns the pattern of the enzyme's digestion
     *
     * @return
     */
    protected abstract Pattern getRestrictionPattern();


    /**
     * Digests the given protein string using the enzyme's restriction pattern
     * and minimal and maximal length allowing no missed cleavages
     *
     * @param protein
     * @return
     */
    public final List<String> digestProtein(String protein, int minLength, int maxLength) {
        return digestProtein(protein, minLength, maxLength, 0);
    }


    /**
     * Digests the given protein string using the enzyme's restriction pattern
     * and minimal and maximal length allowing the given number of missed
     * cleavages.
     *
     * @param protein
     * @return
     */
    public final List<String> digestProtein(String protein, int minLength, int maxLength, int missedCleavages) {
        List<String> peptideList = new ArrayList<String>(
                Arrays.asList(getRestrictionPattern().split(protein)));

        if ((minLength > 0) || (maxLength > 0) || (missedCleavages > 0)) {
            // remove too short or too long peptides and concatenate for missed cleavages
            ListIterator<String> listIt = peptideList.listIterator();
            boolean pepOk;
            List<String> lastPeptides = new ArrayList<String>(missedCleavages);

            while (listIt.hasNext()) {
                String pep = listIt.next();
                pepOk = true;
                if (minLength > 0) {
                    pepOk = (pep.length() >= minLength);
                }

                if (pepOk && (maxLength > 0)) {
                    pepOk = (pep.length() <= maxLength);
                }

                if (!pepOk) {
                    listIt.remove();
                }

                if (missedCleavages > 0) {
                    if (lastPeptides.size() > 0) {
                        StringBuilder uncleavedPep = new StringBuilder(pep);
                        uncleavedPep.reverse();

                        ListIterator<String> lastIt = lastPeptides.listIterator(lastPeptides.size());
                        while (lastIt.hasPrevious()) {
                            uncleavedPep.append(
                                    new StringBuilder(lastIt.previous()).reverse());

                            uncleavedPep.reverse();
                            pepOk = true;
                            if (minLength > 0) {
                                pepOk = (uncleavedPep.length() >= minLength);
                            }
                            if (pepOk && (maxLength > 0)) {
                                pepOk = (uncleavedPep.length() <= maxLength);
                            }
                            if (pepOk) {
                                listIt.add(uncleavedPep.toString());
                            }

                            uncleavedPep.reverse();
                        }
                    }

                    while (lastPeptides.size() >= missedCleavages) {
                        lastPeptides.remove(0);
                    }

                    lastPeptides.add(pep);
                }
            }
        }

        return peptideList;
    }
}