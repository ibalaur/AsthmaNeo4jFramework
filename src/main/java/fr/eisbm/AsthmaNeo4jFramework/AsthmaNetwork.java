package fr.eisbm.AsthmaNeo4jFramework;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;

import fr.eisbm.AsthmaNeo4jFramework.App.LabelTypes;
import fr.eisbm.AsthmaNeo4jFramework.App.RelTypes;

public class AsthmaNetwork {

	private final String ASTHMA_SEEDS_FILE = "files/asthma_seeds.txt";
	private final String HUMAN_TISSUE_ATLAS = "files/HumanTissueAtlas.txt";
	private final String UNIPROT_REACTOME_HOMOSAPIENS_MAPPING = "files/UniProt2Reactome_HomoSapiens.txt";
	private final String GENE_DISEASE_ASSOC_DISGENET_FILE = "files/curated_gene_disease_associations_DisGeNet.txt";
	private final String GENE_UNIPROT_ID_ASSOC_DISGENET_FILE = "files/gene_uniprotId_mapping_DisGeNet.txt";
	private final String PPI_INTACT_FILE = "files/intact_output_protein_interactions.txt";
	private final String ALL_GENE_SYMBOL_NAME_FILE = "files/all_asthma_names_symbols.tab";
	private final String ALL_TRANSCRIPTS_FILE = "files/all_transcripts.tab";
	private final String UNIPROT_ENSEMBLE_MAPPING = "files/Uniprots/mapping_uniprot_ensemble_cell_types.txt";

	private final String DRUG_SPECIFIC_RELATIONSHIPS_FILE = "files/asthma_seeds_drug_relationships/all_drugs_final_data.txt";
	private final String SEQ_SIM_FILE = "files/Asthma terablast_output 19-1-2018.tab";

	private final String ALVEOLAR_EPITHELIAL_INTERACTOME = "files/interactomes/alveolar_epithelial_cell.tsv";
	private final String BLOOD_VESSEL_ENDOTHELIAL_CELL_INTERACTOME = "files/interactomes/blood_vessel_endothelial_cell.tsv";
	private final String DENDRITIC_CELL_INTERACTOME = "files/interactomes/dendritic_cell_-_monocyte_immature_derived.tsv";
	private final String BRONCHIAL_EPITHELIAL_INTERACTOME = "files/interactomes/bronchial_epithelial_cell.tsv";
	private final String BRONCHIAL_SMOOTH_MUSCLE_INTERACTOME = "files/interactomes/bronchial_smooth_muscle_cell.tsv";
	private final String FIBROBLAST_OF_LYMPHATIC_VESSEL_INTERACTOME = "files/interactomes/fibroblast_of_lymphatic_vessel.tsv";
	private final String FIBROBLAST_OF_TUNICA_ADVENTICA_OF_ARTERY_INTERACTOME = "files/interactomes/fibroblast_of_tunica_adventitia_of_artery.tsv";
	private final String IMMATURE_LANGERHANS_INTERACTOME = "files/interactomes/immature_langerhans_cell.tsv";
	private final String LYMPHOCYTE_OF_B_LINEAGE_INTERACTOME = "files/interactomes/lymphocyte_of_b_lineage.tsv";
	private final String MACROPHAGE_INTERACTOME = "files/interactomes/macrophage.tsv";
	private final String MAST_CELL_INTERACTOME = "files/interactomes/mast_cell.tsv";
	private final String MONOCYTE_INTERACTOME = "files/interactomes/monocyte.tsv";
	private final String NATURAL_KILLER_CELLS_INTERACTOME = "files/interactomes/natural_killer_cell.tsv";
	private final String NEUTROPHIL_INTERACTOME = "files/interactomes/neutrophil.tsv";
	private final String T_CELL_INTERACTOME = "files/interactomes/t_cell.tsv";
	private final String SMALL_AIRWAYS_INTERACTOME = "files/interactomes/small_airway_epithelial_cell.tsv";

	public enum CellType_t {
		AlveolarEpithelial, BloodVesel, BronchialEpithelial, BronchialSmoothMuscle, Dendritic_Cell, Fibroblast_of_Lymphatic_Vessel, Fibroblast_of_Tunica_Adventitia_of_Artery, Immature_Lamgerhans, Lymphocyte_of_B_Lineage, Macrophage, MastCell, Monocyte, Natural_Killer_Cells, Neutrophils, TCell, SmallAirways
	}

	Map<String, Node> asthmaSeedSet = new HashMap<String, Node>();
	Map<String, Node> neighboursSet = new HashMap<String, Node>();
	Map<StringPair, Relationship> vPPIMap = new HashMap<StringPair, Relationship>();
	Map<StringPair, Relationship> vCellSpecificPPIMap = new HashMap<StringPair, Relationship>();

	Map<String, Node> diseaseMap = new HashMap<String, Node>();
	Map<String, Node> drugMap = new HashMap<String, Node>();
	Map<String, Node> pathwayMap = new HashMap<String, Node>();
	Map<String, Node> tissueMap = new HashMap<String, Node>();
	Map<String, String> mEnsembleUniprot = new HashMap<String, String>();
	Map<StringPair, Relationship> vDiseaseRelMap = new HashMap<StringPair, Relationship>();
	Map<StringPair, Relationship> vDrugRelMap = new HashMap<StringPair, Relationship>();
	Map<StringPair, Relationship> vTissueRelMap = new HashMap<StringPair, Relationship>();
	Map<StringPair, Relationship> vPathwayMap = new HashMap<StringPair, Relationship>();
	Map<StringPair, Relationship> vSeqSimMap = new HashMap<StringPair, Relationship>();
	String strApostrophe = "\"";
	private String strNoInfo = "NA";

	private void createSeedProteinNode(String szUniprotId, String szProteinName, String szGeneSymbol,
			String szEnsembleTranscript, String szEnsembl) {
		asthmaSeedSet.put(szUniprotId, App.getGraphInstance().createNode(LabelTypes.SeedProtein));
		asthmaSeedSet.get(szUniprotId).setProperty("Label", LabelTypes.SeedProtein.toString());
		asthmaSeedSet.get(szUniprotId).setProperty("ProteinId", strNoInfo);
		asthmaSeedSet.get(szUniprotId).setProperty("UniprotId", szUniprotId);
		asthmaSeedSet.get(szUniprotId).setProperty("ProteinName", szProteinName);
		setEnsemblId(asthmaSeedSet.get(szUniprotId), szEnsembl);
		setGeneSymbols(asthmaSeedSet.get(szUniprotId), szGeneSymbol);

		asthmaSeedSet.get(szUniprotId).setProperty("WELL_REPLICATED", false);
		asthmaSeedSet.get(szUniprotId).setProperty("BIOMARKER", false);
		asthmaSeedSet.get(szUniprotId).setProperty("GWAS", false);
		asthmaSeedSet.get(szUniprotId).setProperty("KANEKO", false);
		asthmaSeedSet.get(szUniprotId).setProperty("LINKAGE_POSITIONAL_CLONING", false);
	}

	private void createNeighbourProtein(String szUniprotId, String szProteinName, String szGeneSymbol,
			String szEnsembl) {
		neighboursSet.put(szUniprotId, App.getGraphInstance().createNode(LabelTypes.NeighbourProtein));
		neighboursSet.get(szUniprotId).setProperty("Label", LabelTypes.NeighbourProtein.toString());
		neighboursSet.get(szUniprotId).setProperty("UniprotId", szUniprotId);
		neighboursSet.get(szUniprotId).setProperty("ProteinName", szProteinName);
		setEnsemblId(neighboursSet.get(szUniprotId), szEnsembl);
		setGeneSymbols(neighboursSet.get(szUniprotId), szGeneSymbol);
	}

	private void createDiseaseNode(String szDiseaseName) {
		diseaseMap.put(szDiseaseName, App.getGraphInstance().createNode(LabelTypes.Disease));
		diseaseMap.get(szDiseaseName).setProperty("Label", LabelTypes.Disease.toString());
		diseaseMap.get(szDiseaseName).setProperty("DiseaseName", szDiseaseName);
	}

	private void createPathwayNode(String szPathwayId, String szPathwayName) {
		pathwayMap.put(szPathwayId, App.getGraphInstance().createNode(LabelTypes.Pathway));
		pathwayMap.get(szPathwayId).setProperty("Label", LabelTypes.Pathway.toString());
		pathwayMap.get(szPathwayId).setProperty("PathwayName", szPathwayName);
	}

	private void createDrugNode(String szDrugName) {
		drugMap.put(szDrugName, App.getGraphInstance().createNode(LabelTypes.Drug));
		drugMap.get(szDrugName).setProperty("DrugName", szDrugName);
		drugMap.get(szDrugName).setProperty("Label", LabelTypes.Drug.toString());
	}

	private void createTissueNode(String szTissueName) {
		tissueMap.put(szTissueName, App.getGraphInstance().createNode(LabelTypes.Tissue));
		tissueMap.get(szTissueName).setProperty("Label", LabelTypes.Tissue.toString());
		tissueMap.get(szTissueName).setProperty("TissueName", szTissueName);
	}

	public void readCellTypeInteractome() throws IOException {

		try (BufferedReader br = new BufferedReader(new FileReader(UNIPROT_ENSEMBLE_MAPPING))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				String szUniprotId = tokens[1];
				String szEnsemblId = tokens[0];

				if (!mEnsembleUniprot.containsKey(szEnsemblId)) {
					mEnsembleUniprot.put(szEnsemblId, szUniprotId);
				}

				if (asthmaSeedSet.containsKey(szUniprotId)) {
					setEnsemblId(asthmaSeedSet.get(szUniprotId), szEnsemblId);
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		loadCellTypeInfo(ALVEOLAR_EPITHELIAL_INTERACTOME, CellType_t.AlveolarEpithelial);
		loadCellTypeInfo(BLOOD_VESSEL_ENDOTHELIAL_CELL_INTERACTOME, CellType_t.BloodVesel);
		loadCellTypeInfo(BRONCHIAL_EPITHELIAL_INTERACTOME, CellType_t.BronchialEpithelial);
		loadCellTypeInfo(BRONCHIAL_SMOOTH_MUSCLE_INTERACTOME, CellType_t.BronchialSmoothMuscle);
		loadCellTypeInfo(DENDRITIC_CELL_INTERACTOME, CellType_t.Dendritic_Cell);
		loadCellTypeInfo(FIBROBLAST_OF_LYMPHATIC_VESSEL_INTERACTOME, CellType_t.Fibroblast_of_Lymphatic_Vessel);
		loadCellTypeInfo(FIBROBLAST_OF_TUNICA_ADVENTICA_OF_ARTERY_INTERACTOME,
				CellType_t.Fibroblast_of_Tunica_Adventitia_of_Artery);
		loadCellTypeInfo(IMMATURE_LANGERHANS_INTERACTOME, CellType_t.Immature_Lamgerhans);
		loadCellTypeInfo(LYMPHOCYTE_OF_B_LINEAGE_INTERACTOME, CellType_t.Lymphocyte_of_B_Lineage);
		loadCellTypeInfo(MACROPHAGE_INTERACTOME, CellType_t.Macrophage);
		loadCellTypeInfo(MAST_CELL_INTERACTOME, CellType_t.MastCell);
		loadCellTypeInfo(MONOCYTE_INTERACTOME, CellType_t.Monocyte);
		loadCellTypeInfo(NATURAL_KILLER_CELLS_INTERACTOME, CellType_t.Natural_Killer_Cells);
		loadCellTypeInfo(NEUTROPHIL_INTERACTOME, CellType_t.Neutrophils);
		loadCellTypeInfo(T_CELL_INTERACTOME, CellType_t.TCell);
		loadCellTypeInfo(SMALL_AIRWAYS_INTERACTOME, CellType_t.SmallAirways);

	}

	public void loadCellTypeInfo(String szCellTypeInteractomeFile, CellType_t eCellType) throws IOException {
		// System.out.println(eCellType);

		if (mEnsembleUniprot.size() > 0) {
			try (BufferedReader br = new BufferedReader(new FileReader(szCellTypeInteractomeFile))) {
				for (String line; (line = br.readLine()) != null;) {
					// process the line.
					String delims = "[\t]";
					String[] tokens = line.split(delims);

					for (int i = 0; i < tokens.length; i++) {
						tokens[i] = tokens[i].toUpperCase().trim();
					}

					String szFirstEnsemble = tokens[0];
					String szSecondEnsembl = tokens[1];
					double iWeight = Math.round(Double.parseDouble(tokens[2]) * 100.0) / 100.0;

					if (mEnsembleUniprot.containsKey(szFirstEnsemble)
							&& mEnsembleUniprot.containsKey(szSecondEnsembl)) {
						String szFirstUniProt = mEnsembleUniprot.get(szFirstEnsemble);
						String szSecondUniProt = mEnsembleUniprot.get(szSecondEnsembl);

						if ((asthmaSeedSet.containsKey(szFirstUniProt))
								|| (asthmaSeedSet.containsKey(szSecondUniProt))) {

							if ((!asthmaSeedSet.containsKey(szFirstUniProt)
									&& (!neighboursSet.containsKey(szFirstUniProt)))) {
								createNeighbourProtein(szFirstUniProt, strNoInfo, strNoInfo, szFirstEnsemble);
							}

							if ((!asthmaSeedSet.containsKey(szSecondUniProt)
									&& (!neighboursSet.containsKey(szSecondUniProt)))) {
								createNeighbourProtein(szSecondUniProt, strNoInfo, strNoInfo, szSecondEnsembl);
							}

							Node _firstProtein = getProteinNode(szFirstUniProt);
							Node _secondProtein = getProteinNode(szSecondUniProt);

							StringPair strPair = new StringPair(szFirstUniProt, szSecondUniProt);
							StringPair strReversePair = new StringPair(szSecondUniProt, szFirstUniProt);

							if ((null != _firstProtein) && (null != _secondProtein)) {

								if ((!vCellSpecificPPIMap.containsKey(strPair))
										&& (!vCellSpecificPPIMap.containsKey(strReversePair))) {

									vCellSpecificPPIMap.put(strPair, _firstProtein.createRelationshipTo(_secondProtein,
											RelTypes.CELL_SPECIFIC_PPI));
								}

								if (vCellSpecificPPIMap.containsKey(strPair)) {
									vCellSpecificPPIMap.get(strPair).setProperty(eCellType.toString(), iWeight);
								}

								if (vCellSpecificPPIMap.containsKey(strReversePair)) {
									vCellSpecificPPIMap.get(strReversePair).setProperty(eCellType.toString(), iWeight);
								}

							}
						}

					}
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}

	private Node getProteinNode(String szUniProt) {
		if (asthmaSeedSet.containsKey(szUniProt)) {
			return asthmaSeedSet.get(szUniProt);
		} else if (neighboursSet.containsKey(szUniProt)) {
			return neighboursSet.get(szUniProt);
		}
		return null;
	}

	public void readAsthmaSeedGenes() throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(ASTHMA_SEEDS_FILE))) {
			for (String line; (line = br.readLine()) != null;) {

				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				if (tokens.length > 2) {

					String szUniprotId = tokens[0];
					String szGeneSymbol = tokens[1];
					String szProvenience = tokens[2];
					String szProteinName = strNoInfo;

					if (tokens.length == 4) {
						szProteinName = tokens[3];
					}

					if (!asthmaSeedSet.containsKey(szUniprotId)) {
						createSeedProteinNode(szUniprotId, szProteinName, szGeneSymbol, strNoInfo, strNoInfo);
					}

					asthmaSeedSet.get(szUniprotId).setProperty(szProvenience, true);
				}
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public void readDataFromIntact() throws IOException {
		createSeedsConnections(PPI_INTACT_FILE);
		createNeighboursConnections(PPI_INTACT_FILE);
	}

	private void createSeedsConnections(String szFileName) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(szFileName))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				String szFirstProtein = tokens[0];
				String szSecondProtein = tokens[1];
				String szInteractDetectionMethod = tokens[2];
				String szReference = tokens[3];
				String szPublicationId = tokens[4];
				String szRelType = tokens[5];
				String szConfidenceLevelInfo = tokens[6];

				if (!szFirstProtein.equals(szSecondProtein)) {

					if (asthmaSeedSet.containsKey(szFirstProtein)) {
						if (asthmaSeedSet.containsKey(szSecondProtein)) {
							createPPIEdge(asthmaSeedSet.get(szFirstProtein), asthmaSeedSet.get(szSecondProtein),
									szRelType, szInteractDetectionMethod, szReference, szPublicationId,
									szConfidenceLevelInfo);
						} else {
							if (!neighboursSet.containsKey(szSecondProtein)) {
								createNeighbourProtein(szSecondProtein, strNoInfo, strNoInfo, strNoInfo);
							}

							createPPIEdge(asthmaSeedSet.get(szFirstProtein), neighboursSet.get(szSecondProtein),
									szRelType, szInteractDetectionMethod, szReference, szPublicationId,
									szConfidenceLevelInfo);
						}
					} else if (asthmaSeedSet.containsKey(szSecondProtein)) {

						if (!neighboursSet.containsKey(szFirstProtein)) {
							createNeighbourProtein(szFirstProtein, strNoInfo, strNoInfo, strNoInfo);
						}

						createPPIEdge(asthmaSeedSet.get(szSecondProtein), neighboursSet.get(szFirstProtein), szRelType,
								szInteractDetectionMethod, szReference, szPublicationId, szConfidenceLevelInfo);

					}
				}

			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	private void createNeighboursConnections(String szFileName) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(szFileName))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				String szFirstProtein = tokens[0];
				String szSecondProtein = tokens[1];
				String szRelType = tokens[5];
				String szInteractDetectionMethod = tokens[2];
				String szReference = tokens[3];
				String szPublicationId = tokens[4];
				String szConfidenceLevelInfo = tokens[6];

				if (!szFirstProtein.equals(szSecondProtein)) {
					if (neighboursSet.containsKey(szFirstProtein) && neighboursSet.containsKey(szSecondProtein)) {
						createPPIEdge(neighboursSet.get(szFirstProtein), neighboursSet.get(szSecondProtein), szRelType,
								szInteractDetectionMethod, szReference, szPublicationId, szConfidenceLevelInfo);
					}
				}

			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private void createPPIEdge(Node szFirstProtein, Node szSecondProtein, String szRelType,
			String szInteractDetectionMethod, String szReference, String szPublicationId,
			String szConfidenceLevelInfo) {
		StringPair strPair = new StringPair(szFirstProtein.getProperty("UniprotId").toString(),
				szSecondProtein.getProperty("UniprotId").toString());
		StringPair strReversePair = new StringPair(szSecondProtein.getProperty("UniprotId").toString(),
				szFirstProtein.getProperty("UniprotId").toString());

		if ((!vPPIMap.containsKey(strPair)) && (!vPPIMap.containsKey(strReversePair))) {

			vPPIMap.put(strPair,
					szFirstProtein.createRelationshipTo(szSecondProtein, Utils.convertStringToRelType(szRelType)));
			vPPIMap.get(strPair).setProperty("Interaction_method_detection", szInteractDetectionMethod);
			vPPIMap.get(strPair).setProperty("Reference", szReference);
			vPPIMap.get(strPair).setProperty("Publication_identifier", szPublicationId);

			int len = szConfidenceLevelInfo.length();
			int istart = 0;
			if (len >= 4) {
				istart = len - 4;
			} else {
				istart = 0;
			}
			String strToken = szConfidenceLevelInfo.substring(istart, len);

			if (strToken.contains("\"")) {
				strToken = strToken.replace("\"", "");
			}

			if (strToken.startsWith(".")) {
				strToken = "0" + strToken;
			}

			vPPIMap.get(strPair).setProperty("Confidence_level", Double.parseDouble(strToken));

		}
	}

	public void readPathwayInfo() throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(UNIPROT_REACTOME_HOMOSAPIENS_MAPPING))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				boolean bFoundPathway = pathwayMap.containsKey(tokens[1]);

				String semicolonDelim = ";";
				String[] szUniprotIds = tokens[0].split(semicolonDelim);

				for (int i = 0; i < szUniprotIds.length; i++) {
					StringPair strPair = new StringPair(szUniprotIds[i], tokens[1]);

					if (!vPathwayMap.containsKey(strPair)) {
						if (asthmaSeedSet.containsKey(szUniprotIds[i])) {

							if (!bFoundPathway) {
								createPathwayNode(tokens[1], tokens[3]);
							}

							vPathwayMap.put(strPair, asthmaSeedSet.get(szUniprotIds[i])
									.createRelationshipTo(pathwayMap.get(tokens[1]), RelTypes.IN_PATHWAY));
							vPathwayMap.get(strPair).setProperty("Reference", "Reactome");

						}
					}
				}
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public void addSeqSimilarityInfo() throws IOException {

		addSeqSimToSeedProteins();
		addSeqSimToNeighProteins();
	}

	private void addSeqSimToSeedProteins() throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(SEQ_SIM_FILE))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				int firstOccurence = tokens[1].indexOf("|");
				int secondOccurence = tokens[1].indexOf("|", firstOccurence + 1);

				String szFirstProtein = tokens[1].substring(firstOccurence + 1, secondOccurence);

				firstOccurence = tokens[4].indexOf("|");
				secondOccurence = tokens[4].indexOf("|", firstOccurence + 1);
				String szSecondProtein = tokens[4].substring(firstOccurence + 1, secondOccurence);

				if (!szFirstProtein.isEmpty() && !szSecondProtein.isEmpty()) {
					if (!szFirstProtein.equals(szSecondProtein)) {
						double fSimScore = Math.round(Double.parseDouble(tokens[7]) * 100.0) / 100.0;
						String szSimSignificance = tokens[8];
						double fAlignmentLength = Math.round(Double.parseDouble(tokens[9]) * 100.0) / 100.0;

						if (asthmaSeedSet.containsKey(szFirstProtein)) {
							if (asthmaSeedSet.containsKey(szSecondProtein)) {
								addSeqSimAssociation(szFirstProtein, szSecondProtein, fSimScore, szSimSignificance,
										fAlignmentLength);
							} else {
								if (!neighboursSet.containsKey(szSecondProtein)) {
									createNeighbourProtein(szSecondProtein, strNoInfo, strApostrophe, strNoInfo);
								}
								addSeqSimAssociation(szFirstProtein, szSecondProtein, fSimScore, szSimSignificance,
										fAlignmentLength);
							}
						} else {
							if (asthmaSeedSet.containsKey(szSecondProtein)) {
								if (!neighboursSet.containsKey(szFirstProtein)) {
									createNeighbourProtein(szFirstProtein, strNoInfo, strApostrophe, strNoInfo);
								}
								addSeqSimAssociation(szSecondProtein, szFirstProtein, fSimScore, szSimSignificance,
										fAlignmentLength);
							}
						}

					}
				}
			}
		}

		catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	private void addSeqSimAssociation(String szFirstProtein, String szSecondProtein, double fSimScore,
			String szSimSignificance, double fAlignLth) {
		Node _firstProtein = getProteinNode(szFirstProtein);
		Node _secondProtein = getProteinNode(szSecondProtein);

		if ((null != _firstProtein) && (null != _secondProtein)) {
			StringPair strPair = new StringPair(szFirstProtein, szSecondProtein);
			StringPair strReversePair = new StringPair(szSecondProtein, szFirstProtein);

			if ((!vSeqSimMap.containsKey(strPair)) && (!vSeqSimMap.containsKey(strReversePair))) {
				vSeqSimMap.put(strPair, _firstProtein.createRelationshipTo(_secondProtein, RelTypes.SEQ_SIM));

				vSeqSimMap.get(strPair).setProperty("Similarity_Score", fSimScore);
				vSeqSimMap.get(strPair).setProperty("Similarity_Significance", szSimSignificance);
				vSeqSimMap.get(strPair).setProperty("Alignment_Length", fAlignLth);
			}
		}
	}

	private void addSeqSimToNeighProteins() throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(SEQ_SIM_FILE))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}
				int firstOccurence = tokens[1].indexOf("|");
				int secondOccurence = tokens[1].indexOf("|", firstOccurence + 1);

				String szFirstProtein = tokens[1].substring(firstOccurence + 1, secondOccurence);

				firstOccurence = tokens[4].indexOf("|");
				secondOccurence = tokens[4].indexOf("|", firstOccurence + 1);
				String szSecondProtein = tokens[4].substring(firstOccurence + 1, secondOccurence);

				if (!szFirstProtein.isEmpty() && !szSecondProtein.isEmpty()) {
					if (!szFirstProtein.equals(szSecondProtein)) {
						if (neighboursSet.containsKey(szFirstProtein) && neighboursSet.containsKey(szSecondProtein)) {
							double fSimScore = Math.round(Double.parseDouble(tokens[7]) * 100.0) / 100.0;
							String szSimSignificance = tokens[8];
							double fAlignmentLength = Math.round(Double.parseDouble(tokens[9]) * 100.0) / 100.0;

							addSeqSimAssociation(szFirstProtein, szSecondProtein, fSimScore, szSimSignificance,
									fAlignmentLength);
						}
					}
				}
			}
		}

		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public void readDataFromDisGeNETFile() throws IOException {
		loadGeneUniprotIdAssociationInfo(GENE_UNIPROT_ID_ASSOC_DISGENET_FILE);
		readGeneDiseaseAssociationInfo(GENE_DISEASE_ASSOC_DISGENET_FILE);
	}

	private void readGeneDiseaseAssociationInfo(String szInputFileName) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(szInputFileName))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				String szProteinId = tokens[0];
				String szGeneSymbol = tokens[1];
				String szUniprotId = null;
				String szScoreValue = tokens[4];
				String szDiseaseName = tokens[3];
				String szReference = tokens[7];

				for (Entry<String, Node> eProtein : asthmaSeedSet.entrySet()) {
					if (eProtein.getValue().getProperty("ProteinId").toString().equals(szProteinId)) {
						szUniprotId = eProtein.getKey();
						break;
					}
				}

				if (null != szUniprotId) {
					asthmaSeedSet.get(szUniprotId).setProperty("GeneSymbol", szGeneSymbol);
					if (!diseaseMap.containsKey(szDiseaseName)) {
						createDiseaseNode(szDiseaseName);
					}

					StringPair strPair = new StringPair(szUniprotId, szDiseaseName);

					if (!vDiseaseRelMap.containsKey(strPair)) {
						vDiseaseRelMap.put(strPair, asthmaSeedSet.get(szUniprotId).createRelationshipTo(
								diseaseMap.get(szDiseaseName), RelTypes.CURATED_DISEASE_GENE_ASSOCIATION));
						vDiseaseRelMap.get(strPair).setProperty("Disease_confidence_level",
								Double.parseDouble(szScoreValue));
						vDiseaseRelMap.get(strPair).setProperty("Disease_reference", szReference);
					}
				}
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private void loadGeneUniprotIdAssociationInfo(String szInputFileName) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(szInputFileName))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				if (asthmaSeedSet.containsKey(tokens[0])) {
					asthmaSeedSet.get(tokens[0]).setProperty("ProteinId", tokens[1]);
				}
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	public void readHumanTissueAtlasInfo() throws IOException {
		readTissueInfo(HUMAN_TISSUE_ATLAS);
		readUniprotIds_Ensembl(ALL_TRANSCRIPTS_FILE);

		createEnsemblTissueAssociation(HUMAN_TISSUE_ATLAS);
	}

	private void createEnsemblTissueAssociation(String szInputFileName) throws IOException {

		try (BufferedReader br = new BufferedReader(new FileReader(szInputFileName))) {

			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String newLine = line.replaceAll(strApostrophe, "");
				String[] tokens = newLine.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				String[] szTissueTokens = tokens[2].split("[:;]");

				for (Entry<String, Node> eProtein : asthmaSeedSet.entrySet()) {
					if (eProtein.getValue().getProperty("EnsemblId").equals(tokens[1])) {

						for (int i = 0; i < szTissueTokens.length; i = i + 2) {
							for (Entry<String, Node> eTissue : tissueMap.entrySet()) {
								if (eTissue.getKey().equals(szTissueTokens[i])) {
									StringPair strPair = new StringPair(tokens[1], eTissue.getKey());
									if (!vTissueRelMap.containsKey(strPair)) {
										vTissueRelMap.put(strPair, eProtein.getValue()
												.createRelationshipTo(eTissue.getValue(), RelTypes.TISSUE_ENHANCED));
										vTissueRelMap.get(strPair).setProperty("RNA_TS_FPKM_value",
												szTissueTokens[i + 1]);
									}
									break;
								}
							}
						}
						break;
					}
				}
			}
		}
	}

	private void readTissueInfo(String szInputFileName) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(szInputFileName))) {

			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String newLine = line.replaceAll(strApostrophe, "");
				String[] tokens = newLine.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				String[] szTissueTokens = tokens[2].split("[:;]");

				for (int i = 0; i < szTissueTokens.length; i = i + 2) {
					if (!tissueMap.containsKey(szTissueTokens[i])) {
						createTissueNode(szTissueTokens[i]);
					}
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private void readUniprotIds_Ensembl(String szInputFileName) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(szInputFileName))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				String semicolonDelim = ";";
				String[] szUniprotIds = tokens[0].split(semicolonDelim);

				for (int i = 0; i < szUniprotIds.length; i++) {
					String szUniprot = szUniprotIds[i];
					if (asthmaSeedSet.containsKey(szUniprot)) {
						setEnsemblId(asthmaSeedSet.get(szUniprot), tokens[1]);
					}
				}
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public void loadProteinAdditionalInfo() throws IOException {
		lodProteinNameSymbol(ALL_GENE_SYMBOL_NAME_FILE);
		loadTranscripts(ALL_TRANSCRIPTS_FILE);
	}

	private void loadTranscripts(String szInputFileName) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(szInputFileName))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				if (asthmaSeedSet.containsKey(tokens[0])) {
					setEnsemblId(asthmaSeedSet.get(tokens[0]), tokens[1]);
				} else if (neighboursSet.containsKey(tokens[0])) {
					setEnsemblId(neighboursSet.get(tokens[0]), tokens[1]);
				}

			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	private void lodProteinNameSymbol(String szInputFileName) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(szInputFileName))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.
				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				if (asthmaSeedSet.containsKey(tokens[0])) {
					setProteinName(asthmaSeedSet.get(tokens[0]), tokens[1]);
					setGeneSymbols(asthmaSeedSet.get(tokens[0]), tokens[2]);
				} else if (neighboursSet.containsKey(tokens[0])) {
					setProteinName(neighboursSet.get(tokens[0]), tokens[1]);
					setGeneSymbols(neighboursSet.get(tokens[0]), tokens[2]);
				}

			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	private void setProteinName(Node node, String szProteinName) {
		node.setProperty("ProteinName", szProteinName);
	}

	private void setGeneSymbols(Node node, String szGeneSymbols) {
		if (!szGeneSymbols.equals("\"")) {
			node.setProperty("GeneSymbol", szGeneSymbols);
		} else {
			node.setProperty("GeneSymbol", strNoInfo);
		}
	}

	private void setEnsemblId(Node node, String szEnsemblId) {
		if (!szEnsemblId.equals("\"")) {
			node.setProperty("EnsemblId", szEnsemblId);
		} else {
			node.setProperty("EnsemblId", strNoInfo);
		}
	}

	public void readDrugSpecificRelationships() throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(DRUG_SPECIFIC_RELATIONSHIPS_FILE))) {
			for (String line; (line = br.readLine()) != null;) {
				// process the line.

				String delims = "[\t]";
				String[] tokens = line.split(delims);

				for (int i = 0; i < tokens.length; i++) {
					tokens[i] = tokens[i].toUpperCase().trim();
				}

				String szUniprotId = tokens[0];
				String szDrugName = tokens[2];
				String szInteractionType = tokens[3];
				String szSource = tokens[4];
				String szVersion = tokens[7];
				String szPubMedIds = strNoInfo;

				if (tokens.length > 8) {
					szPubMedIds = tokens[8];
				}

				boolean bFoundDrug = drugMap.containsKey(szDrugName);
				StringPair strPair = new StringPair(szUniprotId, szDrugName);
				boolean bFoundRel = vDrugRelMap.containsKey(strPair);

				if (!bFoundRel) {
					if (!bFoundDrug) {
						createDrugNode(szDrugName);
					}

					vDrugRelMap.put(strPair, asthmaSeedSet.get(szUniprotId).createRelationshipTo(
							drugMap.get(szDrugName), Utils.convertDrugRelType(szInteractionType)));

					vDrugRelMap.get(strPair).setProperty("Source", szSource);
					vDrugRelMap.get(strPair).setProperty("Version", szVersion);
					vDrugRelMap.get(strPair).setProperty("PubMedIds", szPubMedIds);
				}

			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}
}
