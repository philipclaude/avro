"use strict";

//
// callback when the user interface is to be initialized (called by wv/render.js)
//
function wvInitUI()
{
    wv.clipped = [];
    wv.terminalON = false;

    // zero out the buttons
    init_controls();

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
    canvas.addEventListener('mousedown',  getMouseDown,    false);
    canvas.addEventListener('mousemove',  getMouseMove,    false);
    canvas.addEventListener('mouseup',    getMouseUp,      false);
    canvas.addEventListener("wheel",      getMouseRoll,    false);
    canvas.addEventListener('mouseout',   mouseLeftCanvas, false);
}

function wvClientConnected(text)
{
  // ask the server for the directory listing
  wv.socketUt.send("getdir");
}

//
// callback when a (text) message is received from the server (called by wv/socket.js)
//
function wvServerMessage(text)
{
  // convert the server message to JSON
  var msg = JSON.parse(text.substring(0,text.length-1));

  // figure out what to do
  var commands = msg["commands"];

  for (var i=0;i<commands.length;i++)
  {
    var data = msg[commands[i]+"-data"];

    if (commands[i]=="make_tree")
    {
      make_tree(data["primitives"]);
    }
    else if (commands[i]=="ls")
    {
      wv.dir = data;
    }
    else if (commands[i]=="colorbar")
    {
      var clo = document.getElementById('colorbar-text-lo');
      var chi = document.getElementById('colorbar-text-hi');

      clo.innerHTML = data[0];
      chi.innerHTML = data[1];
    }
  }
}

function init_controls()
{
  document.getElementById('clip-offset').value = 0.;
  document.getElementById('clip-alpha').value = 0.;
  document.getElementById('clip-beta').value = 0.;
  document.getElementById('flip-clip').checked = false;
  document.getElementById('plot-clip').checked = false;
  document.getElementById('real-time').checked = false;
}

function reset_clip()
{
  document.getElementById('clip-offset').value = 0.;
  document.getElementById('clip-alpha').value = 0.;
  document.getElementById('clip-beta').value = 0.;
  var state;
  if (document.getElementById('plot-clip').checked)
    state = "on";
  else
    state = "off";
  wv.socketUt.send("resetclip|"+state);
}


function get_clip_message()
{
  var d = document.getElementById('clip-offset').value;
  var alpha = document.getElementById('clip-alpha').value;
  var beta = document.getElementById('clip-beta').value;
  return String(d)+","+String(alpha)+","+String(beta);
}

function request_clip()
{
  // get all the plot indexes which are to be clipped
  var msg = "clip|";

  // append the clipping plane information
  msg += get_clip_message();

  // send the message to the server
  wv.socketUt.send(msg);
}

function plot_clip()
{
  if (document.getElementById('real-time').checked)
  {
    // also inform the plotter to do the clipping
    request_clip();
    return;
  }

  if (document.getElementById('plot-clip').checked)
  {
    var dir;
    if (document.getElementById('flip-clip').checked)
      dir = "R";
    else dir = "L";
    var msg = "plotclip|"+String(dir)+"|";
    msg += get_clip_message();

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
  print_message("toggle field!!")
  wv.socketUt.send("modgprim|"+iplot+"|"+rank);
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
  var node  = 1;
  var count = 0;
  var plots = [];
  for (var i=0;i<primitives.length;i++)
  {
    var dummy = "";

    var root = node;
    plots.push(root);

    // create a root node
    myTree.addNode( 0 , "Plot "+i , "",dummy,null,"","viz","on",toggle_viz,"edg","on",toggle_edges,"alpha","off",toggle_transparency );
    node++;

    // add the volumes
    var vol_node = node;
    myTree.addNode( root , "Volumes" , "",dummy,null,"","viz","on",toggle_viz,"edg","on",toggle_edges,"alpha","off",toggle_transparency );
    node++;
    var volumes = JSON.parse(primitives[i])["Volumes"];
    for (var j=0;j<volumes.length;j++)
    {
      var gprim = volumes[j];
      myTree.addNode( vol_node , "Volume " + j.toString() , "" , gprim , null,"", "viz","on",toggle_viz,"edg","on",toggle_edges,"alpha","off",toggle_transparency);
      node++;
    }

    // add the faces
    var face_node = node;
    myTree.addNode( root , "Faces" , "",dummy,null,"","viz","on",toggle_viz,"edg","on",toggle_edges,"alpha","off",toggle_transparency );
    node++;
    var faces = JSON.parse(primitives[i])["Faces"];
    for (var j=0;j<faces.length;j++)
    {
      var gprim = faces[j];
      myTree.addNode( face_node , "Face " + j.toString() , "" , gprim , null,"", "viz","on",toggle_viz,"edg","on",toggle_edges,"alpha","off",toggle_transparency);
      node++;
    }

    // add the edges
    var edge_node = node;
    myTree.addNode( root , "Edges" , "",dummy,null,"","viz","on",toggle_viz,"edg","on",toggle_edges );
    node++;
    var edges = JSON.parse(primitives[i])["Edges"];
    for (var j=0;j<edges.length;j++)
    {
      var gprim = edges[j];
      myTree.addNode( edge_node , "Edge " + j.toString() , "" , gprim , null,"", "Viz","on",toggle_viz,"Grd","on" );
      node++;
    }

    // add the nodes
    var node_node = node;
    myTree.addNode( root , "Nodes" , "",dummy,null,"","Viz","on",toggle_viz );
    node++;
    var nodes = JSON.parse(primitives[i])["Nodes"];
    for (var j=0;j<nodes.length;j++)
    {
      var gprim = nodes[j];
      myTree.addNode( node_node , "Node " +j.toString() , "" , gprim , null,"", "Viz","on",toggle_viz );
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

function delete_plots()
{
  myTree.clear();
  myTree.build();
}

function clear_plots()
{
  wv.socketUt.send("clear");
  //wv.sceneGraph = {};
  wvInitUI();
}

//
// callback to toggle Viz property
//
function toggle_viz(e)
{
    if (wv.curMode != 0)
    {
      alert("Command disabled.  Press 'Cancel' or 'OK' first");
      return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Viz property
    if (myTree.valu1[inode] == "off")
    {
      myTree.prop(inode, 1, "on");
    }
    else if (myTree.valu1[inode] == "on")
    {
      myTree.prop(inode, 1, "off");
    }
    else
    {
      alert("Illegal Viz property:"+myTree.valu1[inode]);
      return;
    }
}


//
// callback to toggle Grd property
//
function toggle_edges(e)
{
    if (wv.curMode != 0)
    {
      alert("Command disabled.  Press 'Cancel' or 'OK' first");
      return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the edges property
    if (myTree.valu2[inode] == "off")
    {
      myTree.prop(inode, 2, "on");
    }
    else if (myTree.valu2[inode] == "on")
    {
      myTree.prop(inode, 2, "off");
    }
    else
    {
      alert("Illegal Grd property:"+myTree.valu2[inode]);
      return;
    }
}


//
// callback to toggle Trn property
//
function toggle_transparency(e)
{
    if (wv.curMode != 0)
    {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Trn property
    if (myTree.valu3[inode] == "off")
    {
      myTree.prop(inode, 3, "on");
    }
    else if (myTree.valu3[inode] == "on")
    {
      myTree.prop(inode, 3, "off");
    }
    else
    {
      alert("Illegal Trn property:"+myTree.valu3[inode]);
      return;
    }
}


//
// callback when "onresize" event occurs
//
// resize the frames (with special handling to width of tlframe and height of brframe)
//
function resize_frames()
{
    var scrollSize = 24;

    // get the size of the client (minus amount to account for scrollbars)
    var body = document.getElementById("mainBody");
    var bodyWidth  = body.clientWidth  - scrollSize;
    var bodyHeight = body.clientHeight - scrollSize;

    // get the elements associated with the frames and the canvas
    var topleft = document.getElementById("tlframe");
    var butnfrm = document.getElementById("butnfrm");
    var treefrm = document.getElementById("treefrm");
    var toprite = document.getElementById("trframe");
    var botrite = document.getElementById("brframe");
    var canvas  = document.getElementById(wv.canvasID);
    var plot_collapse = document.getElementById("plot_collapse");
    //plot_collapse.style.height = "600px";

    // compute and set the widths of the frames
    //    (do not make tlframe larger than 250px)
    var leftWidth = Math.round(0.2 * bodyWidth);
    if (leftWidth<350) leftWidth = 350;
    var riteWidth = bodyWidth - leftWidth;
    var canvWidth = riteWidth - scrollSize;

    topleft.style.width = leftWidth+"px";
    treefrm.style.width = "100%";//leftWidth-scrollSize+"px";
    toprite.style.width = riteWidth+"px";
    botrite.style.width = bodyWidth+"px";//riteWidth+"px";
    canvas.style.width  = canvWidth+"px";
    canvas.width        = canvWidth;

    // compute and set the heights of the frames
    var botmHeight = Math.round(0.20 * bodyHeight);
    var  topHeight = bodyHeight - botmHeight;
    var canvHeight =  topHeight -scrollSize;/* - scrollSize - 5;*/

    topleft.style.height =  topHeight+"px";
    treefrm.style.height = "100%";//plot_collapse.style.height;
    toprite.style.height =  topHeight+"px";
    botrite.style.height = botmHeight+"px";
    canvas.style.height  = canvHeight+"px";
    canvas.height        = canvHeight;
}

//
// post a message into the brframe
//
function print_message(mesg,col,newline,brackets)
{
    if (col==undefined) col = "yellow";
    var term = document.getElementById("brframe");

    var pre;
    var end;

    if (newline==undefined) newline = true;

    if (brackets || brackets==undefined)
    {
      pre = "<< ";
      end = " >>";
    }
    else
    {
      pre = "";
      end = "";
    }

    term.print_message(pre+mesg+end,col,newline);
}
