"use strict";

//
// callback when the user interface is to be initialized (called by wv-render.js)
//
function wvInitUI()
{
    wv.clipped = [];
    wv.terminalON = false;

    // zero out the buttons
    initControls();

    deletePlots();

    // set up extra storage for matrix-matrix multiplies
    wv.uiMatrix   = new J3DIMatrix4();
    wv.saveMatrix = new J3DIMatrix4(wv.mvMatrix);

                                   // ui cursor variables
    wv.cursorX   = -1;             // current cursor position
    wv.cursorY   = -1;
    wv.keyPress  = -1;             // last key pressed
    wv.keyCode   = -1;
    wv.startX    = -1;             // start of dragging position
    wv.startY    = -1;
    wv.button    = -1;             // button pressed
    wv.modifier  =  0;             // modifier (shift,alt,cntl) bitflag
    wv.flying    =  1;             // flying multiplier (do not set to 0)
    wv.offTop0   =  0;             // offset to upper-left corner of the canvas
    wv.offLeft0  =  0;
    wv.offTop7   =  0;             // offset to upper-left corner of Sketcher
    wv.offLeft7  =  0;
    wv.dragging  =  false;         // true during drag operation
    wv.picking   =  0;             // keycode of command that turned picking on
    wv.locating  =  0;             // keycode of command that turned locating on
    wv.focus     = [0, 0, 0, 1];   // focus data needed in locating
    wv.builtTo   = 99999;          // last Branch in previous successful build
    wv.menuEvent = undefined;      // event associated with click in Tree
    wv.server    = undefined;      // string passed back from "identify;"
    wv.loLimit   = -1;             // lower limit in key
    wv.upLimit   = +1;             // upper limit in key
    wv.curMode   =  0;             //-1 to disable buttons, etc.
                                   // 0 show WebViewer in canvas
    wv.sgData    = {};             // scene graph metadata
    wv.scale     =  1;             // scale factor for axes
    wv.getFocus  = undefined;      // entry to get first focus

    document.addEventListener('keypress', getKeyPress,     false);
    document.addEventListener('keydown',  getKeyDown,      false);
    document.addEventListener('keyup',    getKeyUp,        false);

    var canvas = document.getElementById(wv.canvasID);
    canvas.addEventListener('mousedown',  getMouseDown0,    false);
    canvas.addEventListener('mousemove',  getMouseMove0,    false);
    canvas.addEventListener('mouseup',    getMouseUp0,      false);
    canvas.addEventListener("wheel",      getMouseRoll0,    false);
    canvas.addEventListener('mouseout',   mouseLeftCanvas0, false);

    //var keycan = document.getElementById(wv.canvasKY);
    //keycan.addEventListener('mouseup',    setKeyLimits,    false);

}

function initControls()
{
  document.getElementById('clipd').value = 0.;
  document.getElementById('clipalpha').value = 0.;
  document.getElementById('clipbeta').value = 0.;
  document.getElementById('flipclip').checked = false;
  document.getElementById('plotclip').checked = false;
  document.getElementById('realtime').checked = false;
}

function resetClip()
{
  document.getElementById('clipd').value = 0.;
  document.getElementById('clipalpha').value = 0.;
  document.getElementById('clipbeta').value = 0.;
  var state;
  if (document.getElementById('plotclip').checked)
    state = "on";
  else
    state = "off";
  wv.socketUt.send("resetclip|"+state);
}

function addTreePlot(iparent,gprim,name)
{
  var graphic = wv.sceneGraph[gprim];
  if (graphic==undefined)
  {
    console.log(wv.sceneGraph);
    console.log("unknown gprim "+gprim+" with name "+name);
    return;
  }
  if (graphic.GPtype==2)
  {
    // triangles
    myTree.addNode(iparent,name, "",gprim,null,"","Viz","on",toggleViz,"Grd","on",toggleGrd,"Trn","off",toggleTrn);
  }
  else if (graphic.GPtype==1)
  {
    // lines
    myTree.addNode(iparent,name, "",gprim,null,"","Viz","on",toggleViz,"Grd","on",toggleGrd,"Trn","off",toggleTrn);
  }
  else if (graphic.GPtype==0)
  {
    // points
    myTree.addNode(iparent,name, "",gprim,null,"","Viz","on",toggleViz);
  }
  myTree.build();
}

function wvClientConnected(text)
{
  // ask the server for the directory listing
  wv.socketUt.send("getdir");
}

function clearOptions(select)
{
  var n = select.options.length;
  for (var i=0;i<n;i++)
  {
    // remove 0'th item because once we do, the next one is the 0th
    select.remove(0);
  }
}

function clippingPlaneMessage()
{
  var d = document.getElementById('clipd').value;
  var alpha = document.getElementById('clipalpha').value;
  var beta = document.getElementById('clipbeta').value;
  return String(d)+","+String(alpha)+","+String(beta);
}

function clippedPlots()
{
  var msg="";
  for (var k=0;k<wv.clipped.length;k++)
  {
    if (wv.clipped[k].checked)
    {
      msg += String(k) +",";
    }
  }
  if (document.getElementById('flipclip').checked) msg += "R";
  else msg += "L";
  msg += "|";
  return msg;
}

function clipRequest()
{
  // get all the plot indexes which are to be clipped
  var msg = "clip|";
  msg += clippedPlots();

  // send the clipping plane information
  msg += clippingPlaneMessage();
  wv.socketUt.send(msg);
}

function plotClip()
{
  if (document.getElementById('realtime').checked)
  {
    // also inform the plotter to do the clipping
    clipRequest();
    return;
  }

  if (document.getElementById('plotclip').checked)
  {
    var dir;
    if (document.getElementById('flipclip').checked)
      dir = "R";
    else dir = "L";
    var msg = "plotclip|"+String(dir)+"|";
    msg += clippingPlaneMessage();

    wv.socketUt.send(msg);
  }
  else
    wv.socketUt.send("hideclip|");

}

function toggle_field(e,select,gprim,iplot)
{
  // get the rank of the field we need to plot
  var rank = select.selectedIndex -1; // -1 to account for 'none' option

  // ask the server to modify the graphics primitives
  postMessage("toggle field!!")
  wv.socketUt.send("modgprim|"+iplot+"|"+rank);
}

function makeClipbox( iplot , inode )
{
  // get the html element corresponding to this table row
  var elem = document.getElementById("node"+inode);

  // append a select to the end of this row
  var input = document.createElement("input");
  input.type = "checkbox";
  input.id = "clip"+iplot;
  input.name = "clip";
  input.value = true;
  input.title = "whether this plot is clipped by the clipping plane";

  wv.clipped.push( input );

  elem.appendChild(input);

}



function addPlot()
{
  var option = document.getElementById("directory");
  if (option.selectedIndex==0) return;
  name = option.options[option.selectedIndex].text;

  postMessage("sending request for "+name);

  // we made it here, so the server should send back a mesh
  wv.socketUt.send("cd|"+name);

}

function makePlots(info)
{
  var ptr    = info["ptr"];
  var parent = info["parent"];
  var name   = info["name"];

  // number of table entries
  var nt = ptr.length;

  for (var i=0;i<nt;i++)
    addTreePlot( parent[i] , ptr[i] , name[i] );

  wvUpdateUI();

  // add the dropdowns after the UI has been updated because otherwise the tree drawing overwrites the dropdown
  var iplot = 0;
  for (var i=0;i<nt;i++)
  {
    // we need to add a drop down for the fields if this is made in the root node (i.e. a plot)
    if (parent[i]==0)
    {
      // make a select button for the field on this plot
      makeSelect(iplot,i+2,ptr[i],info); // remember html table is offset by 2

      // make a clip checkbox for this plot
      makeClipbox(iplot,i+2,ptr[i],info);

      iplot++;
    }
  }
}

function addOption(select,text,value)
{
  var opt = document.createElement("option");
  opt.text = text;
  opt.value = value;
  select.options.add(opt);
}

function updateDirectoryList(info)
{
  var entries = info["entries"];
  var drop = document.getElementById("directory");
  clearOptions(drop);
  for (var i=0;i<entries.length;i++)
    addOption(drop,entries[i],"");
}

function addPlotDropdown(info)
{
  var entries = info["entries"];

  // the first entry is the plot number, the second is the name
  var iplot = entries[0];
  var name = entries[1];
  var select = document.getElementById("select"+iplot);
  addOption(select,name,"");
  select.selectedIndex = select.options.length -1;
}

function make_fields( node , fields )
{
  // get the html element corresponding to this table row
  var elem = document.getElementById("node"+node);

  // append a select to the end of this row
  var select = document.createElement("select");
  select.id = "select"+node;

  // first add a 'none' option to only see the geometry
  var option0 = document.createElement("option");
  option0.text = "none";
  option0.value = "";
  select.options.add(option0);

  // add the remaining field options
  var nf = fields.length;
  for (var i=0;i<nf;i++)
  {
    var option1 = document.createElement("option");
    option1.text = fields[i];
    option1.value = "";

    select.options.add(option1);
  }

  if (select.options.length==1)
    select.selectedIndex = 0;
  else
    select.selectedIndex = 1; // the server will colour the field in the first rank

  select.onchange = function(e) { toggle_field(e,select,node); };
  elem.appendChild(select);
}

function make_tree(primitives)
{
  var node  = 2; // offset by two because of tree header
  var count = 0;
  var plots = [];
  for (var i=0;i<primitives.length;i++)
  {
    var dummy = "";

    var root = node;
    plots.push(root);

    // create a root node
    myTree.addNode( 0 , "mesh"+i , "",dummy,null,"","Viz","on",toggleViz,"Grd","on",toggleGrd,"Trn","off",toggleTrn );
    node++;

    // add the volumes
    var vol_node = node;
    myTree.addNode( root , "Volumes" , "",dummy,null,"","Viz","on",toggleViz,"Grd","on",toggleGrd,"Trn","off",toggleTrn );
    node++;
    var volumes = JSON.parse(primitives[i])["Volumes"];
    for (var j=0;j<volumes.length;j++)
    {
      var gprim = volumes[j];
      myTree.addNode( vol_node , gprim , "" , gprim , null,"", "Viz","on",toggleViz,"Grd","on",toggleGrd,"Trn","off",toggleTrn);
      node++;
    }

    // add the faces
    var face_node = node;
    myTree.addNode( root , "Faces" , "",dummy,null,"","Viz","on",toggleViz,"Grd","on",toggleGrd,"Trn","off",toggleTrn );
    node++;
    var faces = JSON.parse(primitives[i])["Faces"];
    for (var j=0;j<faces.length;j++)
    {
      var gprim = faces[j];
      myTree.addNode( face_node , gprim , "" , gprim , null,"", "Viz","on",toggleViz,"Grd","on",toggleGrd,"Trn","off",toggleTrn);
      node++;
    }

    // add the edges
    var edge_node = node;
    myTree.addNode( root , "Edges" , "",dummy,null,"","Viz","on",toggleViz,"Grd","on",toggleGrd );
    node++;
    var edges = JSON.parse(primitives[i])["Edges"];
    for (var j=0;j<edges.length;j++)
    {
      var gprim = edges[j];
      myTree.addNode( edge_node , gprim , "" , gprim , null,"", "Viz","on",toggleViz,"Grd","on" );
      node++;
    }

    // add the nodes
    var node_node = node;
    myTree.addNode( root , "Nodes" , "",dummy,null,"","Viz","on",toggleViz );
    node++;
    var nodes = JSON.parse(primitives[i])["Nodes"];
    for (var j=0;j<nodes.length;j++)
    {
      var gprim = nodes[j];
      myTree.addNode( node_node , gprim , "" , gprim , null,"", "Viz","on",toggleViz );
      node++;
    }

  }

  // build the tree (must be done before adding fields in order to create node id's)
  myTree.build();

  for (var i=0;i<plots.length;i++)
  {
    var fields = JSON.parse(primitives[i])["fields"];
    make_fields( plots[i] , fields );
  }
}

//
// callback when a (text) message is received from the server (called by wv-socket.js)
//
function wvServerMessage(text)
{
  // convert the server message to JSON
  var msg = JSON.parse(text.substring(0,text.length-1));

  // figure out what to do
  var commands = msg["commands"];

  for (var i=0;i<commands.length;i++)
  {
    postMessage(commands[i]);
    var data = msg[commands[i]+"-data"];

    if (commands[i]=="make_tree")
    {
      make_tree(data["primitives"]);
    }
  }

  return;
}

function deletePlots()
{
  myTree.clear();
  myTree.addNode(0, "Plots", "show plots", "", null, "*"); // offsets by 2
  myTree.build();
}

function clearPlots()
{
  wv.socketUt.send("clear");
  //wv.sceneGraph = {};
  wvInitUI();
}

function rollbackDirectory()
{
  wv.socketUt.send("cd|"+"..");
}
