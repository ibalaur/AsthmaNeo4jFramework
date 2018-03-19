package fr.eisbm.AsthmaNeo4jFramework;

import fr.eisbm.AsthmaNeo4jFramework.App.RelTypes;

public class Utils {
	static public RelTypes convertStringToRelType(String string) {

		if ((string.toUpperCase().trim().contains("PHYSICAL ASSOCIATION"))
				|| (string.toUpperCase().trim().contains("ASSOCIATION"))
				|| (string.toUpperCase().trim().contains("PHYSICAL INTERACTION"))
				|| (string.toUpperCase().trim().contains("DIRECT INTERACTION"))
				|| (string.toUpperCase().trim().contains("COVALENT BINDING"))
				|| (string.toUpperCase().trim().contains("ENZYMATIC REACTION"))
				|| (string.toUpperCase().trim().contains("SELF INTERACTION"))
				|| (string.toUpperCase().trim().contains("DISULFIDE BOND"))
				|| (string.toUpperCase().trim().contains("TRANSGLUTAMINATION REACTION"))
				|| (string.toUpperCase().trim().contains("PHOSPHORYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("SUMOYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("CLEAVAGE REACTION"))
				|| (string.toUpperCase().trim().contains("DEPHOSPHORYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("UBIQUITINATION REACTION"))
				|| (string.toUpperCase().trim().contains("OXIDOREDUCTASE ACTIVITY ELECTRON TRANSFER REACTION"))
				|| (string.toUpperCase().trim().contains("ACETYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("DEUBIQUITINATION REACTION"))
				|| (string.toUpperCase().trim().contains("GTPASE REACTION"))
				|| (string.toUpperCase().trim().contains("PROTEIN CLEAVAGE"))
				|| (string.toUpperCase().trim().contains("DEACETYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("METHYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("PALMITOYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("NEDDYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("HYDROXYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("ADP RIBOSYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("PROLINE ISOMERIZATION  REACTION"))
				|| (string.toUpperCase().trim().contains("PHOSPHOTRANSFER REACTION"))
				|| (string.toUpperCase().trim().contains("DEMETHYLATION REACTION"))
				|| (string.toUpperCase().trim().contains("RNA CLEAVAGE"))
				|| (string.toUpperCase().trim().contains("DENEDDYLATION REACTION")))
			return RelTypes.PPI_ASSOCIATION;

		if (string.toUpperCase().trim().contains("COLOCALIZATION"))
			return RelTypes.PPI_COLOCALIZATION;

		if (string.toUpperCase().trim().contains("GENETIC INTERACTION"))
			return RelTypes.PPI_GENETIC_INTERACTION;

		if (string.toUpperCase().trim().contains("PREDICTED INTERACTION"))
			return RelTypes.PPI_PREDICTED_INTERACTION;

		return RelTypes.eNoEvent;
	}

	public static RelTypes convertDrugRelType(String szInteractionType) {
		if (szInteractionType.toUpperCase().trim().equals("ANTAGONIST"))
			return RelTypes.ANTAGONIST;
		if (szInteractionType.toUpperCase().trim().equals("ACTIVATOR"))
			return RelTypes.ACTIVATOR;
		if (szInteractionType.toUpperCase().trim().equals("AGONIST"))
			return RelTypes.AGONIST;
		if (szInteractionType.toUpperCase().trim().equals("ANTIBODY"))
			return RelTypes.ANTIBODY;
		if (szInteractionType.toUpperCase().trim().equals("ANTISENSE"))
			return RelTypes.ANTISENSE;
		if (szInteractionType.toUpperCase().trim().equals("BINDER"))
			return RelTypes.BINDER;
		if (szInteractionType.toUpperCase().trim().equals("INDUCER"))
			return RelTypes.INDUCER;
		if (szInteractionType.toUpperCase().trim().equals("INHIBITOR"))
			return RelTypes.INHIBITOR;
		if (szInteractionType.toUpperCase().trim().equals("MODULATOR"))
			return RelTypes.MODULATOR;
		if (szInteractionType.toUpperCase().trim().equals("POTENTIATOR"))
			return RelTypes.POTENTIATOR;
		if (szInteractionType.toUpperCase().trim().equals("STIMULATOR"))
			return RelTypes.STIMULATOR;
		if (szInteractionType.toUpperCase().trim().equals("SUPPRESSOR"))
			return RelTypes.SUPPRESSOR;
		if (szInteractionType.toUpperCase().trim().equals("ANTAGONIST/PARTIAL AGONIST"))
			return RelTypes.ANTAGONIST_PARTIAL_AGONIST;
		if (szInteractionType.toUpperCase().trim().equals("NEGATIVE MODULATOR"))
			return RelTypes.NEGATIVE_MODULATOR;
		if (szInteractionType.toUpperCase().trim().equals("PARTIAL AGONIST"))
			return RelTypes.PARTIAL_AGONIST;

		return RelTypes.eNoEvent;
	}
}
