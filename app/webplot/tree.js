//
// constructor for a Tree
//
function Tree(doc, treeId)
{
    // alert("in Tree(doc="+doc+", treeId="+treeId+")");

    // remember the document
    this.document = doc;
    this.treeId   = treeId;

    // arrays to hold the Nodes
    this.name    = new Array();
    this.tooltip = new Array();
    this.gprim   = new Array();
    this.click   = new Array();
    this.cmenu   = new Array();
    this.parent  = new Array();
    this.child   = new Array();
    this.next    = new Array();
    this.nprop   = new Array();
    this.opened  = new Array();

    this.prop1  = new Array();
    this.valu1  = new Array();
    this.cbck1  = new Array();
    this.prop2  = new Array();
    this.valu2  = new Array();
    this.cbck2  = new Array();
    this.prop3  = new Array();
    this.valu3  = new Array();
    this.cbck3  = new Array();

    // initialize Node=0 (the root)
    this.name[  0]  = "**root**";
    this.tooltip[0] = "";
    this.gprim[ 0]  = "";
    this.click[ 0]  = null;
    this.cmenu[ 0]  = "";
    this.parent[0]  = -1;
    this.child[ 0]  = -1;
    this.next[  0]  = -1;
    this.nprop[ 0]  =  0;
    this.prop1[ 0]  = "";
    this.valu1[ 0]  = "";
    this.cbck1[ 0]  = null;
    this.prop2[ 0]  = "";
    this.valu2[ 0]  = "";
    this.cbck2[ 0]  = null;
    this.prop3[ 0]  = "";
    this.valu3[ 0]  = "";
    this.cbck3[ 0]  = null;
    this.opened[0]  = +1;

    // add methods
    this.addNode  = TreeAddNode;
    this.expand   = TreeExpand;
    this.contract = TreeContract;
    this.prop     = TreeProp;
    this.clear    = TreeClear;
    this.build    = TreeBuild;
    this.update   = TreeUpdate;
}


//
// add a Node to the Tree
//
function TreeAddNode(iparent, name, tooltip, gprim, click, cmenu, prop1, valu1, cbck1,
                     prop2, valu2, cbck2, prop3, valu3, cbck3)
{
    // alert("in TreeAddNode(iparent="+iparent+", name="+name+", tooltip="+tooltip+", gprim="+gprim+", click="+click+
    //                       ", cmenu="+cmenu+", prop1="+prop1+", valu1="+valu1+",cbck1="+cbck1+
    //                       ", prop2="+prop2+", valu2="+valu2+", cbck2="+cbck2+", prop3="+prop3+
    //                       ", valu3="+valu3+", cbck3="+cbck3+")");

    // validate the input
    if (iparent < 0 || iparent >= this.name.length) {
        alert("iparent="+iparent+" is out of range");
        return;
    }

    // find the next Node index
    var inode = this.name.length;

    // store this Node's values
    this.name[   inode] = name;
    this.tooltip[inode] = tooltip;
    this.gprim[  inode] = gprim;
    this.click[  inode] = click;
    this.cmenu[  inode] = cmenu;
    this.parent[ inode] = iparent;
    this.child[  inode] = -1;
    this.next[   inode] = -1;
    this.nprop[  inode] =  0;
    this.opened[ inode] =  0;

    // store the properties
    if (prop1 !== undefined) {
        this.nprop[inode] = 1;
        this.prop1[inode] = prop1;
        this.valu1[inode] = valu1;
        this.cbck1[inode] = cbck1;
    }

    if (prop2 !== undefined) {
        this.nprop[inode] = 2;
        this.prop2[inode] = prop2;
        this.valu2[inode] = valu2;
        this.cbck2[inode] = cbck2;
    }

    if (prop3 !== undefined) {
        this.nprop[inode] = 3;
        this.prop3[inode] = prop3;
        this.valu3[inode] = valu3;
        this.cbck3[inode] = cbck3;
    }

    // if the parent does not have a child, link this
    //    new Node to the parent
    if (this.child[iparent] < 0) {
        this.child[iparent] = inode;

    // otherwise link this Node to the last parent's child
    } else {
        var jnode = this.child[iparent];
        while (this.next[jnode] >= 0) {
            jnode = this.next[jnode];
        }

        this.next[jnode] = inode;
    }
}


//
// build the Tree (ie, create the html table from the Nodes)
//
function TreeBuild()
{
    var doc = this.document;

    // if the table already exists, delete it and all its children (3 levels)
    var thisTable = doc.getElementById(this.treeId);
    if (thisTable) {
        var child1 = thisTable.lastChild;
        while (child1) {
            var child2 = child1.lastChild;
            while (child2) {
                var child3 = child2.lastChild;
                while (child3) {
                    child2.removeChild(child3);
                    child3 = child2.lastChild;
                }
                child1.removeChild(child2);
                child2 = child1.lastChild;
            }
            thisTable.removeChild(child1);
            child1 = thisTable.lastChild;
        }
        thisTable.parentNode.removeChild(thisTable);
    }

    // build the new table
    var newTable = doc.createElement("table");
    newTable.setAttribute("id", this.treeId);
    doc.getElementById("treefrm").appendChild(newTable);

    newTable.style.borderColor = "white";
    newTable.style.width = doc.getElementById("treefrm").style.width;
    newTable.style.borderRight = "0px";
    doc.getElementById("treefrm").style.borderRight = "0px";

    // traverse the Nodes using depth-first search
    var inode = 1;
    while (inode > 0) {

        // table row "node"+inode
        var newTR = doc.createElement("TR");
        newTR.setAttribute("id", "node"+inode);
        newTable.appendChild(newTR);
        newTable.borderColor = "white";

        // table data "node"+inode+"col1"
        var newTDcol1 = doc.createElement("TD");
        newTDcol1.setAttribute("id", "node"+inode+"col1");
        newTR.appendChild(newTDcol1);

        var newTexta = doc.createTextNode("");
        newTDcol1.appendChild(newTexta);

        // table data "node"+inode+"col2"
        var newTDcol2 = doc.createElement("TD");
        newTDcol2.setAttribute("id", "node"+inode+"col2");
        if (this.cmenu[inode] != "") {
            newTDcol2.className = "fakelinkcmenu";
            if (this.tooltip[inode].length > 0) {
                newTDcol2.title = this.tooltip[inode];
            }
        }
        newTR.appendChild(newTDcol2);

        var newTextb = doc.createTextNode(this.name[inode]);
        newTDcol2.appendChild(newTextb);

        var name = this.name[inode].replace(/\u00a0/g, "");

        // table data "node"+inode+"col3"
        if (this.nprop[inode] > 0) {
            var newTDcol3 = doc.createElement("TD");
            newTDcol3.setAttribute("id", "node"+inode+"col3");
            if (this.cbck1[inode] != "") {
                newTDcol3.className = "fakelinkon";
            }
            newTR.appendChild(newTDcol3);

            if (this.nprop[inode] == 1) {
                newTDcol3.setAttribute("colspan", "3");
            }

            var newTextc = doc.createTextNode(this.prop1[inode]);
            newTDcol3.appendChild(newTextc);
        }

        // table data "node:+inode+"col4"
        if (this.nprop[inode] > 1) {
            var newTDcol4 = doc.createElement("TD");
            newTDcol4.setAttribute("id", "node"+inode+"col4");
            if (this.cbck2[inode] != "") {
                newTDcol4.className = "fakelinkon";
            }
            newTR.appendChild(newTDcol4);

            if (this.nprop[inode] == 2) {
                newTDcol4.setAttribute("colspan", "2");
            }

            var newTextd = doc.createTextNode(this.prop2[inode]);
            newTDcol4.appendChild(newTextd);
        }

        // table data "node:+inode+"col5"
        if (this.nprop[inode] > 2) {
            var newTDcol5 = doc.createElement("TD");
            newTDcol5.setAttribute("id", "node"+inode+"col5");
            if (this.cbck3[inode] != "") {
                newTDcol5.className = "fakelinkon";
            }
            newTR.appendChild(newTDcol5);

            var newTextd = doc.createTextNode(this.prop3[inode]);
            newTDcol5.appendChild(newTextd);
        }

        // go to next row
        if        (this.child[inode] >= 0) {
            inode = this.child[inode];
        } else if (this.next[inode] >= 0) {
            inode = this.next[inode];
        } else {
            while (inode > 0) {
                inode = this.parent[inode];
                if (this.parent[inode] == 0) {
                    newTR = doc.createElement("TR");
                    newTR.setAttribute("height", "10px");
                    newTable.appendChild(newTR);
                }
                if (this.next[inode] >= 0) {
                    inode = this.next[inode];
                    break;
                }
            }
        }
    }

    this.update();
    resize_frames();
}


//
// clear the Tree
//
function TreeClear()
{
    // remove all but the first Node
    this.name.splice(   1);
    this.tooltip.splice(1);
    this.gprim.splice(  1);
    this.click.splice(  1);
    this.cmenu.splice(  1);
    this.parent.splice( 1);
    this.child.splice(  1);
    this.next.splice(   1);
    this.nprop.splice(  1);
    this.opened.splice( 1);

    this.prop1.splice(1);
    this.valu1.splice(1);
    this.cbck1.splice(1);
    this.prop2.splice(1);
    this.valu2.splice(1);
    this.cbck2.splice(1);
    this.prop3.splice(1);
    this.valu3.splice(1);
    this.cbck3.splice(1);

    // reset the root Node
    this.parent[0] = -1;
    this.child[ 0] = -1;
    this.next[  0] = -1;
}


//
// expand a Node in the Tree
//
function TreeContract(inode)
{
    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    }

    // contract inode
    this.opened[inode] = 0;

    // contract all descendents of inode
    for (var jnode = 1; jnode < this.parent.length; jnode++) {
        var iparent = this.parent[jnode];
        while (iparent > 0) {
            if (iparent == inode) {
                this.opened[jnode] = 0;
                break;
            }

            iparent = this.parent[iparent];
        }
    }

    // update the display
    this.update();
}


//
// expand a Node in the Tree
//
function TreeExpand(inode)
{
    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    }

    // expand inode
    this.opened[inode] = 1;

    // update the display
    this.update();
}


//
// change a property of a Node
//
function TreeProp(inode, iprop, onoff)
{
    // alert("in TreeProp(inode="+inode+", iprop="+iprop+", onoff="+onoff+")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    } else if (onoff != "on" && onoff != "off") {
        alert("onoff="+onoff+" is not 'on' or 'off'");
        return;
    }

    // set the property for inode
    if (iprop == 1) {
        this.valu1[inode] = onoff;
    } else if (iprop == 2) {
        this.valu2[inode] = onoff;
    } else if (iprop == 3) {
        this.valu3[inode] = onoff;
    } else {
        alert("iprop="+iprop+" is not 1, 2, or 3");
        return;
    }

    // set property of all descendents of inode
    for (var jnode = 1; jnode < this.parent.length; jnode++) {
        var iparent = this.parent[jnode];
        while (iparent > 0) {
            if (iparent == inode) {
                if        (iprop == 1) {
                    this.valu1[jnode] = onoff;
                } else if (iprop == 2) {
                    this.valu2[jnode] = onoff;
                } else if (iprop == 3) {
                    this.valu3[jnode] = onoff;
                }
                break;
            }

            iparent = this.parent[iparent];
        }
    }

    this.update();
}


//
// update the Tree (after build/expension/contraction/property-set)
//
function TreeUpdate()
{
    var doc = this.document;

    // traverse the Nodes using depth-first search
    for (var inode = 1; inode < this.opened.length; inode++) {
        var element = doc.getElementById("node"+inode);

        // unhide the row
        element.style.display = "table-row";

        // hide the row if one of its parents has .opened=0
        var jnode = this.parent[inode];
        while (jnode != 0) {
            if (this.opened[jnode] == 0) {
                element.style.display = "none";
                break;
            }

            jnode = this.parent[jnode];
        }

        // if the current Node has children, set up appropriate event handler to expand/collapse
        if (this.child[inode] > 0) {
            if (this.opened[inode] == 0) {
                var myElem = doc.getElementById("node"+inode+"col1");
                var This   = this;

                myElem.className = "fakelinkexpand";
                myElem.firstChild.nodeValue = "+";
                myElem.title   = "Expand";
                myElem.onclick = function () {
                    var thisNode = this.id.substring(4);
                    thisNode     = thisNode.substring(0,thisNode.length-4);
                    This.expand(thisNode);
                };

            } else {
                var myElem = doc.getElementById("node"+inode+"col1");
                var This   = this;

                myElem.className = "fakelinkexpand";
                myElem.firstChild.nodeValue = "-";
                myElem.title   = "Collapse";
                myElem.onclick = function () {
                    var thisNode = this.id.substring(4);
                    thisNode     = thisNode.substring(0,thisNode.length-4);
                    This.contract(thisNode);
                };
            }
        }

        if (this.click[inode] !== null) {
            var myElem = doc.getElementById("node"+inode+"col2");
            myElem.onclick = this.click[inode];
        }

        // set the class of the properties
        if (this.nprop[inode] >= 1) {
            var myElem = doc.getElementById("node"+inode+"col3");
            myElem.onclick = this.cbck1[inode];

            if (this.prop1[inode] == "viz") {
                if (this.valu1[inode] == "off") {
                    myElem.setAttribute("class",   "fakelinkoff");
                    myElem.title = "Toggle visibility";
                    if (this.gprim[inode] != "") {
                      if (wv.sceneGraph[this.gprim[inode]]!=undefined)
                      {
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.ON;
                        wv.sceneUpd = 1;
                      }
                    }
                } else {
                    myElem.setAttribute("class",   "fakelinkon");
                    myElem.title = "Toggle visibility";
                    if (this.gprim[inode] != "") {
                        //this.gprim[inode] = "tess";

                        if (wv.sceneGraph[this.gprim[inode]]!=undefined)
                        {
                          wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.ON;
                          wv.sceneUpd = 1;
                        }
                    }
                }
            }
        }

        if (this.nprop[inode] >= 2) {
            var myElem = doc.getElementById("node"+inode+"col4");
            myElem.onclick = this.cbck2[inode];

            if (this.prop2[inode] == "edg") {
                if (this.valu2[ inode] == "off") {
                    myElem.setAttribute("class",   "fakelinkoff");
                    myElem.title = "Toggle edges";

                    if (this.gprim[inode] != "") {
                      if (wv.sceneGraph[this.gprim[inode]]!=undefined)
                      {
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.LINES;
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.POINTS;
                        wv.sceneUpd = 1;
                      }
                    }
                } else {
                    myElem.setAttribute("class",   "fakelinkon");
                    myElem.title = "Toggle edges";

                    if (this.gprim[inode] != "") {
                      if (wv.sceneGraph[this.gprim[inode]]!=undefined)
                      {
                        wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.LINES;
                        wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.POINTS;
                        wv.sceneUpd = 1;
                      }
                    }
                }
            }
        }

        if (this.nprop[inode] >= 3) {
            var myElem = doc.getElementById("node"+inode+"col5");
            myElem.onclick = this.cbck3[inode];

            if (this.prop3[inode] == "alpha") {
                if (this.valu3[ inode] == "off") {
                    myElem.setAttribute("class",   "fakelinkoff");
                    myElem.title = "Toggle transparency";

                    if (this.gprim[inode] != "") {
                      if (wv.sceneGraph[this.gprim[inode]]!=undefined)
                      {
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.TRANSPARENT;
                        wv.sceneUpd = 1;
                      }
                    }
                } else {
                    myElem.setAttribute("class",   "fakelinkon");
                    myElem.title = "Toggle transparency";

                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.TRANSPARENT;
                        wv.sceneUpd = 1;
                    }
                }
            }
        }
    }
}
