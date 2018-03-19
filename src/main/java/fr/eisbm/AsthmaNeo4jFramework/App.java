package fr.eisbm.AsthmaNeo4jFramework;

import java.io.File;
import java.io.IOException;

import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.io.fs.FileUtils;

public class App 
{
	// START SNIPPET: vars
		private static final File databaseDirectory = new File("target/graph.db");

		private static long start = 0;
		private static long end = 0;
		private static GraphDatabaseService graphDb;

		private AsthmaNetwork asthmaGraph = new AsthmaNetwork();

		// END SNIPPET: vars

		// START SNIPPET: createReltype
		public static enum RelTypes implements RelationshipType {
			eNoEvent, IN_PATHWAY, SEQ_SIM, CURATED_DISEASE_GENE_ASSOCIATION, PPI_ASSOCIATION, PPI_COLOCALIZATION, PPI_GENETIC_INTERACTION, PPI_PREDICTED_INTERACTION, TISSUE_ENHANCED, CELL_SPECIFIC_PPI, ANTAGONIST, ACTIVATOR, AGONIST, ANTAGONIST_PARTIAL_AGONIST, ANTIBODY, ANTISENSE, BINDER, INDUCER, INHIBITOR, MODULATOR, NEGATIVE_MODULATOR, PARTIAL_AGONIST, POTENTIATOR, STIMULATOR, SUPPRESSOR
		}

		public static enum LabelTypes implements Label {
			Pathway, Disease, SeedProtein, NeighbourProtein, Tissue, Drug
		}

		public static GraphDatabaseService getGraphInstance() {
			return graphDb;
		}

		// END SNIPPET: createReltype

		public static void main(final String[] args) {
			App asthmaGraphDB = new App();
			try {
				asthmaGraphDB.createDb();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			asthmaGraphDB.runEmbeddedQueries_SharedProperties_normalised(true);

			shutDown();
		}

		private void runEmbeddedQueries_SharedProperties_normalised(boolean b_SharedProteins) {

			String query_create_shared_disease = "MATCH (p1:SeedProtein)--(d:Disease)--(p2:SeedProtein) WHERE id(p1) > id(p2) "
					+ " WITH p1,p2, COUNT(d) AS comm "
					+ " OPTIONAL MATCH (p1)--(d1:Disease) WITH p1,p2, comm, COUNT(d1) AS x1 "
					+ " OPTIONAL MATCH (p2)--(d2:Disease) WITH p1,p2, comm, x1, COUNT(d2) AS x2 "
					+ " UNWIND [x1,x2] AS val WITH p1,p2,comm, MAX(val) AS maximum"
					+ " CREATE (p1)-[:SHARED_DISEASES{disease_weight: ROUND(100 * comm/maximum) / 100 }]->(p2)";

			createSharedProperties(query_create_shared_disease);

			String query_create_shared_pathway = "MATCH (p1:SeedProtein)--(path:Pathway)--(p2:SeedProtein) WHERE id(p1) > id(p2) "
					+ " WITH p1,p2, COUNT(path) AS comm "
					+ " OPTIONAL MATCH (p1)--(path1:Pathway) WITH p1,p2, comm, COUNT(path1) AS x1 "
					+ " OPTIONAL MATCH (p2)--(path2:Pathway) WITH p1,p2, comm, x1, COUNT(path2) AS x2 "
					+ " UNWIND [x1,x2] AS val WITH p1,p2,comm, MAX(val) AS maximum"
					+ " CREATE (p1)-[:SHARED_PATHWAYS{pathway_weight: ROUND(100 * comm/maximum) / 100 }]->(p2)";

			createSharedProperties(query_create_shared_pathway);

			String query_create_shared_tissue = "MATCH (p1:SeedProtein)--(t:Tissue)--(p2:SeedProtein) WHERE id(p1) > id(p2) "
					+ " WITH p1,p2, COUNT(t) AS comm "
					+ " OPTIONAL MATCH (p1)--(t1:Tissue) WITH p1,p2, comm, COUNT(t1) AS x1 "
					+ " OPTIONAL MATCH (p2)--(t2:Tissue) WITH p1,p2, comm, x1, COUNT(t2) AS x2 "
					+ " UNWIND [x1,x2] AS val WITH p1,p2,comm, MAX(val) AS maximum"
					+ " CREATE (p1)-[:SHARED_TISSUES{tissue_weight: ROUND(100 * comm/maximum) / 100 }]->(p2)";

			createSharedProperties(query_create_shared_tissue);

			String query_create_shared_drug = "MATCH (p1:SeedProtein)--(drug:Drug)--(p2:SeedProtein) WHERE id(p1) > id(p2) "
					+ " WITH p1,p2, COUNT(drug) AS comm "
					+ " OPTIONAL MATCH (p1)--(drug1:Drug) WITH p1,p2, comm, COUNT(drug1) AS x1 "
					+ " OPTIONAL MATCH (p2)--(drug2:Drug) WITH p1,p2, comm, x1, COUNT(drug2) AS x2 "
					+ " UNWIND [x1,x2] AS val WITH p1,p2,comm, MAX(val) AS maximum"
					+ " CREATE (p1)-[:SHARED_DRUGS{drug_weight: ROUND(100 * comm/maximum) / 100 }]->(p2)";

			createSharedProperties(query_create_shared_drug);

			if (true == b_SharedProteins) {

				String query_create_shared_prot = "MATCH (p1:SeedProtein)--(neigh:NeighbourProtein)--(p2:SeedProtein) WHERE id(p1) > id(p2) "
						+ " WITH p1,p2, COUNT(neigh) AS comm "
						+ " OPTIONAL MATCH (p1)--(neigh1:NeighbourProtein) WITH p1,p2, comm, COUNT(neigh1) AS x1 "
						+ " OPTIONAL MATCH (p2)--(neigh2:NeighbourProtein) WITH p1,p2, comm, x1, COUNT(neigh2) AS x2 "
						+ " UNWIND [x1,x2] AS val WITH p1,p2,comm, MAX(val) AS maximum"
						+ " CREATE (p1)-[:SHARED_PROTEINS{protein_weight: ROUND(100 * comm/maximum) / 100 }]->(p2)";

				createSharedProperties(query_create_shared_prot);
			}
		}

		private void createSharedProperties(String query) {

			try (Transaction tx = graphDb.beginTx(); Result result = graphDb.execute(query)) {
				tx.success();
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}

		void createDb() throws IOException {

			FileUtils.deleteRecursively(databaseDirectory);

			// START SNIPPET: startDb
			graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(databaseDirectory);
			registerShutdownHook(graphDb);
			// END SNIPPET: startDb

			// START SNIPPET: transaction
			try (Transaction tx = graphDb.beginTx()) {
				// Database operations go here
				// END SNIPPET: transaction

				// START SNIPPET: addData
				try {

					start = System.currentTimeMillis();
					System.out.println("readAsthmaSeedGenes");
					asthmaGraph.readAsthmaSeedGenes();
					System.out.println("readHumanTissueAtlasInfo");
					asthmaGraph.readHumanTissueAtlasInfo();
					System.out.println("loadProteinAdditionalInfo");
					asthmaGraph.loadProteinAdditionalInfo();
					System.out.println("readCellTypeInfo");
					asthmaGraph.readCellTypeInteractome();
					System.out.println("readDataFromIntact");
					asthmaGraph.readDataFromIntact();
					System.out.println("readPathwayInfo");
					asthmaGraph.readPathwayInfo();
					System.out.println("readDrugSpecificRelationships");
					asthmaGraph.readDrugSpecificRelationships();
					System.out.println("readDataFromDisGeNETFile");
					asthmaGraph.readDataFromDisGeNETFile();
					System.out.println("addSeqSimilarityInfo");
					asthmaGraph.addSeqSimilarityInfo();

					end = System.currentTimeMillis();
					long diff = end - start;
					System.out.println("Time (in seconds) is : " + diff * 0.001);

				} catch (IOException e) {
					e.printStackTrace();
				} catch (Exception e) {
					e.printStackTrace();
				}

				tx.success();
			}
			// END SNIPPET: transaction

		}

		static void shutDown() {
			System.out.println();
			System.out.println("Shutting down database ...");

			// START SNIPPET: shutdownServer

			try {
				graphDb.shutdown();
			} catch (Exception e) {
				e.printStackTrace();
			}
			// END SNIPPET: shutdownServer
		}

		// START SNIPPET: shutdownHook

		private static void registerShutdownHook(final GraphDatabaseService graphDb) {

			// Registers a shutdown hook for the Neo4j instance so that it shuts
			// down nicely when the VM exits (even if you "Ctrl-C" the

			// running application).

			Runtime.getRuntime().addShutdownHook(new Thread() {
				@Override
				public void run() {
					graphDb.shutdown();
				}
			});

		}

		// END SNIPPET: shutdownHook
		private static void deleteFileOrDirectory(File file) {
			if (file.exists()) {
				if (file.isDirectory()) {
					for (File child : file.listFiles()) {
						deleteFileOrDirectory(child);
					}
				}
				file.delete();
			}
		}
}
