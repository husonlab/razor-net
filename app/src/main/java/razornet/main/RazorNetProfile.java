/*
 * RazorNetProfile.java Copyright (C) 2026 Daniel H. Huson
 *
 *  (Some files contain contributions from other authors, who are then mentioned separately.)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

package razornet.main;

import splitstree6.algorithms.network.network2view.ShowNetwork;
import splitstree6.algorithms.source.source2characters.CharactersLoader;
import splitstree6.algorithms.source.source2distances.DistancesLoader;
import splitstree6.algorithms.source.source2network.NetworkLoader;
import splitstree6.data.*;
import splitstree6.main.IAppProfile;
import splitstree6.workflow.Algorithm;
import splitstree6.workflow.AlgorithmNode;
import splitstree6.workflow.DataBlock;
import splitstree6.workflow.Workflow;

import java.util.Set;
import java.util.function.Predicate;

/**
 * The application profile of the RazorNet program: the algorithms that it offers, the default workflow for each
 * type of input, and the menu items that it keeps.
 */
public class RazorNetProfile implements IAppProfile {
	private static final Set<String> WHITELIST = Set.of(
			// Loaders
			"CharactersLoader", "DistancesLoader", "NetworkLoader", "TaxaFilter",
			// Filters
			"CharactersTaxaFilter", "DistancesTaxaFilter", "CharactersFilter",
			"NetworkTaxaFilter", "StretchFilter", "SimplifyFilter",
			// The RazorNet algorithms
			"RazorNet", "RazorHaploNet", "TCS", "TCSHaplotypeNetwork",
			// related:
			"NeighborNet", "MedianJoining", "MinSpanningNetwork",
			// Viewer
			"ShowNetwork"
	);

	@Override
	public String getName() {
		return "RazorNet";
	}

	@Override
	public Predicate<Algorithm<?, ?>> getAlgorithmFilter() {
		return a -> WHITELIST.contains(a.getClass().getSimpleName());
	}

	@Override
	public boolean setupWorkflow(Workflow workflow, Class<? extends DataBlock> inputType, SourceBlock sourceBlock) {
		if (inputType.equals(CharactersBlock.class)) {
			workflow.setupInputAndWorkingNodes(sourceBlock, new CharactersLoader(),
					new TaxaBlock(), new CharactersBlock());
			var networkNode = workflow.newDataNode(new NetworkBlock());
			workflow.newAlgorithmNode(new razornet.algorithms.RazorHaploNet(),
					workflow.getWorkingTaxaNode(),
					workflow.getWorkingDataNode(), networkNode);
			var viewerNode = workflow.newDataNode(new ViewBlock());
			workflow.newAlgorithmNode(new ShowNetwork(),
					workflow.getWorkingTaxaNode(), networkNode, viewerNode);
			return true;
		}
		if (inputType.equals(DistancesBlock.class)) {
			workflow.setupInputAndWorkingNodes(sourceBlock, new DistancesLoader(),
					new TaxaBlock(), new DistancesBlock());
			var networkNode = workflow.newDataNode(new NetworkBlock());
			workflow.newAlgorithmNode(new razornet.algorithms.RazorNet(),
					workflow.getWorkingTaxaNode(),
					workflow.getWorkingDataNode(), networkNode);
			var viewerNode = workflow.newDataNode(new ViewBlock());
			workflow.newAlgorithmNode(new ShowNetwork(),
					workflow.getWorkingTaxaNode(), networkNode, viewerNode);
			return true;
		}
		if (inputType.equals(NetworkBlock.class)) {
			// a network is opened and shown as it is, without recomputing it
			workflow.setupInputAndWorkingNodes(sourceBlock, new NetworkLoader(),
					new TaxaBlock(), new NetworkBlock());
			var viewerNode = workflow.newDataNode(new ViewBlock());
			workflow.newAlgorithmNode(new ShowNetwork(),
					workflow.getWorkingTaxaNode(), workflow.getWorkingDataNode(), viewerNode);
			return true;
		}
		return false;   // everything else uses SplitsTree6's default
	}

	@Override
	public Predicate<Class<? extends DataBlock>> getOpenableInputTypeFilter() {
		return t -> t == CharactersBlock.class || t == DistancesBlock.class || t == NetworkBlock.class;
	}

	@Override
	public boolean shouldReplaceWorkflow(Workflow loadedWorkflow) {
		// Keep as-is only if the workflow already contains a RazorNet algorithm.
		return loadedWorkflow.nodes().stream()
				.filter(n -> n instanceof AlgorithmNode)
				.map(n -> ((AlgorithmNode<?, ?>) n).getAlgorithm())
				.noneMatch(a -> a instanceof razornet.algorithms.RazorNet || a instanceof razornet.algorithms.RazorHaploNet);
	}

	private static final Set<String> KEEP_MENU_IDS = Set.of("newMenuItem", "openMenuItem", "openRecentMenu", "replaceDataMenuItem",
			"editInputMenuItem", "importDialogMenuItem", "importHaplotypeMenuItem", "importTaxonDisplayMenuItem", "importTaxonTraitsMenuItem",
			"saveMenuItem", "saveAsMenuItem", "exportImageMenuItem", "exportTaxonDisplayLabelsMenuItem", "exportTaxonTraitsMenuItem",
			"exportWorkflowMenuItem", "pageSetupMenuItem", "printMenuItem", "closeMenuItem", "quitMenuItem",
			"undoMenuItem", "redoMenuItem", "cutMenuItem", "copyMenuItem", "copyImageMenuItem",
			"pasteMenuItem", "duplicateMenuItem", "deleteMenuItem", "findAgainMenuItem", "gotoLineMenuItem",
			"selectAllMenuItem", "selectNoneMenuItem", "selectInverseMenuItem", "selectFromPreviousMenuItem", "selectBracketsMenuItem",
			"selectCompatibleSitesMenuItem", "addTaxonSetMenuItem", "removeTaxonSetMenuItem", "useFullScreenMenuItem", "increaseFontSizeMenuItem",
			"decreaseFontSizeMenuItem", "zoomInMenuItem", "zoomOutMenuItem", "zoomInHorizontalMenuItem", "zoomOutHorizontalMenuItem",
			"resetMenuItem", "rotateLeftMenuItem", "rotateRightMenuItem", "flipMenuItem", "layoutLabelsMenuItem",
			"filterTaxaMenuItem", "filterCharactersMenuItem", "GroupIdenticalHaplotypesMenuItem", "traitsMenuItem",
			"setWindowSizeMenuItem", "checkForUpdatesMenuItem", "aboutMenuItem");
	;

	@Override
	public Set<String> getKeepMenuIds() {
		return KEEP_MENU_IDS;
	}
}
