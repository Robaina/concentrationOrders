let nodeSize = 200;
let nodeFontSize = 200;
let edgeFontSize = 150;
let edgeWidth = 10;
let fontColor = '#ba9938';
let childrenColor = '#5a61c2';
let parentColor = '#ed5e9c';
let selectedColor = '#17cfad';
let backgroundColor = 'rgb(54, 54, 54)';

graphStyle = [
  {
    selector: 'node',
    style: {
        'label': 'data(label)',
        'width': nodeSize + '%',
        'height': nodeSize + '%',
        'color': fontColor,
        'background-color': '#9e9e9e',
        'font-size': 0,
        'text-halign': 'center'
    }
  },

  {
    selector: 'edge',
    style: {
        'label': 'data(label)',
        'color': fontColor,
        'font-size': 0,
        'width': edgeWidth + '%',
        'line-color': 'grey',
        'target-arrow-color': 'grey',
        'target-arrow-shape': 'triangle',
        'control-point-step-size': '140px',
        'curve-style': 'unbundled-bezier'
    }
  },

  {
    selector: '.childrenEdges',
    style: {
        'font-size': edgeFontSize,
        'width': '4%',
        'line-color': childrenColor,
        'target-arrow-color': childrenColor,
        'arrow-scale': 3
    }
  },

  {
    selector: '.parentEdges',
    style: {
        'font-size': edgeFontSize,
        'width': '4%',
        'line-color': parentColor,
        'target-arrow-color': parentColor,
        'arrow-scale': 3
    }
  },

  {
    selector: '.childrenNodes',
    style: {
      'background-color': childrenColor,
      'text-background-color': backgroundColor,
      'text-background-opacity': 0,
      'font-size': nodeFontSize,
    }
  },

  {
    selector: '.parentNodes',
    style: {
      'background-color': parentColor,
      'text-background-color': backgroundColor,
      'text-background-opacity': 0,
      'font-size': nodeFontSize,
    }
  },

  {
    selector: '.selectedNode',
    style: {
      'background-color': selectedColor,
      'z-index': 1,
      'font-size': nodeFontSize,
    }
  },

  {
    selector: '.selectedSubGraphNode',
    style: {
      'label': 'data(name)',
      'background-color': selectedColor,
      'text-background-color': backgroundColor,
      'text-background-opacity': 1,
      'z-index': 1,
      'font-size': nodeFontSize,
    }
  }

]
