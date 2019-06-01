/* Flux order relation graph of Escherichia coli
   Semidán Robaina Estévez, 2018
*/
let rxnNameLabels = rxnLabels = [];
let option;
let selectedNodeID = "#nad_c";
let aboutButtonPressed = false;
let button, selectedContainer, buttonPressed;
let oldLabel, oldColor;
let inputRxnName, inputRxnLabel, hoveredNodeID, subcyHoveredNodeID;
let graphContainer = document.getElementById("cy");


for (node of data['nodes']) {

  inputRxnName = node['data']['name'];
  inputRxnLabel = node['data']['label'];
  rxnNameLabels.push(inputRxnName);
  rxnLabels.push(inputRxnLabel);
  option = document.createElement('option');
  option.text = inputRxnName + " (" + inputRxnLabel + ")";
  option.value = inputRxnLabel;
  option.setAttribute("id", inputRxnLabel);
  if ("#" + inputRxnLabel === selectedNodeID) {
    option.selected = "selected";
  }
  document.getElementById("select-list").add(option);
}

function changeSelectedReaction() {
  let selector = document.getElementById("select-list");
  selectedNodeID = "#" + selector[selector.selectedIndex].value;
  initializeGraph(selectedNodeID);
}

// Define graph object
let cy = cytoscape({
  container: graphContainer,
  elements: data,
  style: graphStyle,
  userZoomingEnabled: false,
  autoungrabify: false,
  userPanningEnabled: false
});

// Initialize graph
function initializeGraph(selectedNodeID) {

  selectNodes(cy, selectedNodeID);
  cy.layout(options).run();
  // plotChart(selectedNodeID);

  // Modify the position of some nodes a little bit
//   xposGLYCDx = cy.$('#GLYCDx').renderedPosition('x');
//   yposGLYCDx = cy.$('#GLYCDx').renderedPosition('y');
//   xposGLYK = cy.$('#GLYK').renderedPosition('x');
//   yposGLYK = cy.$('#GLYK').renderedPosition('y');
//   cy.$('#GLYCDx').renderedPosition('x', (1 - 0.5) * xposGLYCDx);
//   cy.$('#GLYCDx').renderedPosition('y', (1 + 0.08) * yposGLYCDx);
//   cy.$('#GLYK').renderedPosition('x', (1 + 0.45) * xposGLYK);
//   cy.$('#GLYK').renderedPosition('y', (1 + 0.08) * yposGLYK);
//   xposATPS4rpp = cy.$('#ATPS4rpp').renderedPosition('x');
//   cy.$('#ATPS4rpp').renderedPosition('x', xposATPS4rpp - 50);
//   xposACONTa = cy.$('#ACONTa').renderedPosition('x');
//   cy.$('#ACONTa').renderedPosition('x', xposACONTa + 25);
//   xposACONTb = cy.$('#ACONTb').renderedPosition('x');
//   cy.$('#ACONTb').renderedPosition('x', xposACONTb + 50);
}

initializeGraph(selectedNodeID);

// Interactive block
cy.on('mouseover', 'node', function(event) {
  hoveredNodeID = '#' + this.id();
  if (hoveredNodeID !== selectedNodeID) {
    cy.$(hoveredNodeID).addClass('selectedNode');
  }
  cy.$(hoveredNodeID).addClass('selectedSubGraphNode');
});
cy.on('mouseout', 'node', function(event) {
  if (hoveredNodeID !== selectedNodeID) {
    cy.$(hoveredNodeID).removeClass('selectedNode');
  }
  cy.$(hoveredNodeID).removeClass('selectedSubGraphNode');
});

cy.on('click tap', 'node', function(event) {
  selectedNodeID = '#' + this.id();
  document.getElementById(this.id()).selected = "selected";
  selectNodes(cy, selectedNodeID);
  // plotChart(selectedNodeID);
});

function selectNodes(cy, selectedNodeID) {
  childrenNodes = cy.$(selectedNodeID).successors('node');
  parentNodes = cy.$(selectedNodeID).predecessors('node');
  childrenEdges = cy.$(selectedNodeID).successors('edge');
  parentEdges = cy.$(selectedNodeID).predecessors('edge');

  // Change style classes on click
  cy.batch(function() {
    cy.elements().not(childrenNodes, parentNodes).classes('node');
    childrenNodes.classes('childrenNodes');
    childrenEdges.classes('childrenEdges');
    parentNodes.classes('parentNodes');
    parentEdges.classes('parentEdges');
    cy.$(selectedNodeID).classes('selectedNode');
  });

}

function showContainer(selectedContainer) {

  let containers = document.getElementsByClassName("container");
  for (let i=0; i<containers.length; i++) {
    containers[i].style.display = "none";
  }
  selectedContainer.style.display = "initial";
}

function isMediaScreen() {
  return window.innerWidth < 768;
}

function appendCreditsToAboutIfMedia() {
  credits = document.getElementById("credits");
  credits.style.display = "intitial";
  document.getElementById("about-container").appendChild(credits);
}

if (isMediaScreen()) {
  appendCreditsToAboutIfMedia();
}

// function showAbout() {
//   aboutButtonPressed = !aboutButtonPressed;
//   let container = document.getElementById("about-container");
//
//   if (window.innerWidth < 768) {
//       document.getElementById("reaction-form").style.display = "none";
//       document.getElementById("main-title").innerHTML = "The flux order relation";
//       showContainer(container);
//
//   } else {
//
//     if (aboutButtonPressed) {
//       container.style.display = "initial";
//       setTimeout(() => container.style.opacity = 1, 100);
//     } else {
//       container.style.opacity = 0;
//       setTimeout(() => container.style.display = "none", 1500);
//     }
//   }
// }

// Media screens buttons

// function showPlot() {
//   let container = document.getElementById("plot-container");
//   document.getElementById("main-title").innerHTML = "Metabolic subsystems of the subgraph";
//   document.getElementById("reaction-form").style.display = "none";
//   showContainer(container);
//   plotChart(); //xtick labels lost if not replotted
// }

function showGraph() {
  let container = document.getElementById("graph-container");
  document.getElementById("main-title").innerHTML = "Ordering of metabolic fluxes in <em>Escherichia coli</em>";
  document.getElementById("reaction-form").style.display = "initial";
  showContainer(container);
}

function appendCreditsToAboutIfMedia() {
  credits = document.getElementById("credits")
  document.getElementById("about-container").appendChild(credits);
}

if (isMediaScreen()) {
  appendCreditsToAboutIfMedia();
}

function isMediaScreen() {
  return window.innerWidth < 768;
}
