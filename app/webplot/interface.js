// this part was mostly written by John Dannenhoffer

//
// callback when the user interface should be updated (called by wv/render.js)
//
function wvUpdateUI()
{
    // special code for delayed-picking mode
    if (wv.picking > 0)
    {

        // if something is picked, post a message
        if (wv.picked !== undefined) {

            // second part of '^' operation
            if (wv.picking == 94) {
                var mesg = "Picked: "+wv.picked.gprim;

                try {
                    var attrs = wv.sgData[wv.picked.gprim];
                    for (var i = 0; i < attrs.length; i+=2) {
                        mesg = mesg + "\n        "+attrs[i]+"= "+attrs[i+1];
                    }
                } catch (x) {
                }
                print_message(mesg);
            }

            wv.picked  = undefined;
            wv.picking = 0;
            wv.pick    = 0;

        // abort picking on mouse motion
        } else if (wv.dragging) {
            print_message("Picking aborted");

            wv.picking = 0;
            wv.pick    = 0;
        }

        wv.keyPress = -1;
        wv.dragging = false;
    }

    // special code for delayed-locating mode
    if (wv.locating > 0)
    {

        // if something is located, post a message
        if (wv.located !== undefined) {

            // second part of '@' operation
            if (wv.locating == 64) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                print_message("Located: x="+xloc.toFixed(4)+", y="+yloc.toFixed(4)+
                                   ", z="+zloc.toFixed(4));
            }

            wv.located  = undefined;
            wv.locating = 0;
            wv.locate   = 0;

        // abort locating on mouse motion
        } else if (wv.dragging) {
            print_message("Locating aborted");

            wv.locating = 0;
            wv.locate   = 0;
        }

        wv.keyPress = -1;
        wv.dragging = false;
    }

    // deal with key presses
    if (wv.keyPress != -1 && wv.curMode == 0 && !wv.terminalON)
    {

        var myKeyPress = String.fromCharCode(wv.keyPress);

        // '?' -- help
        if (myKeyPress == "?") {
            print_message("........................... Viewer Cursor options ...........................\n" +
                        "ctrl-h <Home> - initial view             ctrl-f        - front view          \n" +
                        "ctrl-l        - leftside view            ctrl-r        - riteside view       \n" +
                        "ctrl-t        - top view                 ctrl-b        - bottom view         \n" +
                        "ctrl-i <PgUp> - zoom in                  ctrl-o <PgDn> - zoom out            \n" +
                        "<Left>        - rotate or xlate left     <Rite>        - rotate or xlate rite\n" +
                        "<Up>          - rotate or xlate up       <Down>        - rotate or xlate down\n" +
                        ">             - save view                <             - recall view         \n" +
                        "^             - query object at cursor   @             - get coords. @ cursor\n" +
                        "!             - toggle flying mode       *             - center view @ cursor\n" +
                        ".............................. ? - get help .................................");

        // '^' -- query at cursor
        } else if (myKeyPress == '^') {
            wv.picking  = 94;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // '@' -- locate at cursor
        } else if (myKeyPress == '@') {
            wv.locating = 64;
            wv.locate   = 1;

        // '*' -- center view
        } else if (myKeyPress == '*') {
            wv.centerV = 1;
            print_message("View being centered");

        // '!' -- toggle flying mode
        } else if (myKeyPress == "!") {
            if (wv.flying <= 1) {
                print_message("Turning flying mode ON");
                wv.flying = 10;
            } else {
                print_message("Turning flying mode OFF");
                wv.flying = 1;
            }

        // '>' -- save view
        } else if (myKeyPress == ">") {
            print_message("Saving current view");
            wv.saveMatrix.load(wv.mvMatrix);
            wv.sceneUpd = 1;

        // '<' -- recall view
        } else if (myKeyPress == "<") {
            print_message("Restoring saved view");
            wv.mvMatrix.load(wv.saveMatrix);
            wv.sceneUpd = 1;

        // '<Home>' -- initial view
        } else if (wv.keyPress == 0 && wv.keyCode == 36) {
            wv.mvMatrix.makeIdentity();
            wv.scale    = 1;
            wv.sceneUpd = 1;

        // '<end>' -- not used
        } else if (wv.keyPress == 0 && wv.keyCode == 35) {
            print_message("<End> is not supported.  Use '?' for help");

        // '<PgUp>' -- zoom in
        } else if (wv.keyPress == 0 && wv.keyCode == 33) {
            if (wv.modifier == 0) {
                wv.mvMatrix.scale(2.0, 2.0, 2.0);
                wv.scale *= 2.0;
            } else {
                wv.mvMatrix.scale(1.25, 1.25, 1.25);
                wv.scale *= 1.25;
            }
            wv.sceneUpd = 1;

        // '<PgDn>' -- zoom out
        } else if (wv.keyPress == 0 && wv.keyCode == 34) {
            if (wv.modifier == 0) {
                wv.mvMatrix.scale(0.5, 0.5, 0.5);
                wv.scale *= 0.5;
            } else {
                wv.mvMatrix.scale(0.8, 0.8, 0.8);
                wv.scale *= 0.8;
            }
            wv.sceneUpd = 1;

        // '<Delete>' -- not used
        } else if (wv.keyPress == 0 && wv.keyCode == 46) {
            print_message("<Delete> is not supported.  Use '?' for help");

        // '<Left>' -- rotate or translate object left
        } else if (wv.keyPress == 0 && wv.keyCode == 37) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(-30, 0,1,0);
                } else {
                    wv.mvMatrix.rotate( -5, 0,1,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(-0.5, 0.0, 0.0);
                } else {
                    wv.mvMatrix.translate(-0.1, 0.0, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Right>' -- rotate or translate object right
        } else if (wv.keyPress == 0 && wv.keyCode == 39) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(+30, 0,1,0);
                } else {
                    wv.mvMatrix.rotate( +5, 0,1,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(+0.5, 0.0, 0.0);
                } else {
                    wv.mvMatrix.translate(+0.1, 0.0, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Up>' -- rotate or translate object up
        } else if (wv.keyPress == 0 && wv.keyCode == 38) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(-30, 1,0,0);
                } else {
                    wv.mvMatrix.rotate( -5, 1,0,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(0.0, +0.5, 0.0);
                } else {
                    wv.mvMatrix.translate(0.0, +0.1, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Down>' -- rotate or translate object down
        } else if (wv.keyPress == 0 && wv.keyCode == 40) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(+30, 1,0,0);
                } else {
                    wv.mvMatrix.rotate( +5, 1,0,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(0.0, -0.5, 0.0);
                } else {
                    wv.mvMatrix.translate(0.0, -0.1, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // 'ctrl-h' - initial view (same as <Home>)
        } else if ((wv.keyPress == 104 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   8 && wv.modifier == 4 && wv.keyCode ==  8)   ) {
            cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-i' - zoom in (same as <PgUp> without shift)
        } else if ((wv.keyPress == 105 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   9 && wv.modifier == 4 && wv.keyCode ==  9)   ) {
            cmdIn();
            wv.keyPress = -1;
            return;

        // 'ctrl-o' - zoom out (same as <PgDn> without shift)
        } else if ((wv.keyPress == 111 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  15 && wv.modifier == 4 && wv.keyCode == 15)   ) {
            cmdOut();
            wv.keyPress = -1;
            return;

        // 'ctrl-f' - front view (same as <Home>)
        } else if ((wv.keyPress == 102 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   6 && wv.modifier == 4 && wv.keyCode ==  6)   ) {
            cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-r' - riteside view
        } else if ((wv.keyPress == 114 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  18 && wv.modifier == 4 && wv.keyCode == 18)   ) {
            cmdRite();
            wv.keyPress = -1;
            return;

        // 'ctrl-l' - leftside view
        } else if ((wv.keyPress == 108 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  12 && wv.modifier == 4 && wv.keyCode == 12)   ) {
            cmdLeft();
            wv.keyPress = -1;
            return;

        // 'ctrl-t' - top view
        } else if ((wv.keyPress == 116 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  20 && wv.modifier == 4 && wv.keyCode == 20)   ) {
            cmdTop();
            wv.keyPress = -1;
            return;

        // 'ctrl-b' - bottom view
        } else if ((wv.keyPress ==  98 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   2 && wv.modifier == 4 && wv.keyCode ==  2)   ) {
            cmdBotm();
            wv.keyPress = -1;
            return;

            // NOP
        } else if (wv.keyPress == 0 && wv.modifier == 0) {

        // unknown command
        }
        else if ( wv.keyPress==115 && wv.modifier==0  ) // 's' for step
        {
          wv.socketUt.send("step");
        }
        else if ( wv.keyPress==114 && wv.modifier==0 ) // 'r' for rewind
        {
          wv.socketUt.send("rewind");
        }
        else {
            print_message("Unknown command (keyPress="+wv.keyPress
                        +", modifier="+wv.modifier
                        +", keyCode="+wv.keyCode+").  Use '?' for help");
        }

        wv.keyPress = -1;
    }

    // UI is in screen coordinates (not object)
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();

    // deal with mouse movement
    if (wv.dragging)
    {

        // cntrl is down (rotate)
        if (wv.modifier == 4) {
            var angleX =  (wv.startY - wv.cursorY) / 4.0 / wv.flying;
            var angleY = -(wv.startX - wv.cursorX) / 4.0 / wv.flying;
            if ((angleX != 0.0) || (angleY != 0.0)) {
                wv.mvMatrix.rotate(angleX, 1,0,0);
                wv.mvMatrix.rotate(angleY, 0,1,0);
                wv.sceneUpd = 1;
            }

        // alt-shift is down (rotate)
        } else if (wv.modifier == 3) {
            var angleX =  (wv.startY - wv.cursorY) / 4.0 / wv.flying;
            var angleY = -(wv.startX - wv.cursorX) / 4.0 / wv.flying;
            if ((angleX != 0.0) || (angleY != 0.0)) {
                wv.mvMatrix.rotate(angleX, 1,0,0);
                wv.mvMatrix.rotate(angleY, 0,1,0);
                wv.sceneUpd = 1;
            }

        // alt is down (spin)
        } else if (wv.modifier == 2) {
            var xf = wv.startX - wv.width  / 2;
            var yf = wv.startY - wv.height / 2;

            if ((xf != 0.0) || (yf != 0.0)) {
                var theta1 = Math.atan2(yf, xf);
                xf = wv.cursorX - wv.width  / 2;
                yf = wv.cursorY - wv.height / 2;

                if ((xf != 0.0) || (yf != 0.0)) {
                    var dtheta = Math.atan2(yf, xf) - theta1;
                    if (Math.abs(dtheta) < 1.5708) {
                        var angleZ = 128*(dtheta) / 3.1415926 / wv.flying;
                        wv.mvMatrix.rotate(angleZ, 0,0,1);
                        wv.sceneUpd = 1;
                    }
                }
            }

        // shift is down (zoom)
        } else if (wv.modifier == 1) {
            if (wv.cursorY != wv.startY) {
                var scale = Math.exp((wv.cursorY - wv.startY) / 512.0 / wv.flying);
                wv.mvMatrix.scale(scale, scale, scale);
                wv.scale   *= scale;
                wv.sceneUpd = 1;
            }

        // no modifier (translate)
        } else {
            var transX = (wv.cursorX - wv.startX) / 256.0 / wv.flying;
            var transY = (wv.cursorY - wv.startY) / 256.0 / wv.flying;
            if ((transX != 0.0) || (transY != 0.0)) {
                wv.mvMatrix.translate(transX, transY, 0.0);
                wv.sceneUpd = 1;
            }
        }

        // if not flying, then update the start coordinates
        if (wv.flying <= 1) {
            wv.startX = wv.cursorX;
            wv.startY = wv.cursorY;
        }
    }
}

//
// callback when the server goes down either on abort or normal closure (called by wv-socket.js)
//
function wvServerDown()
{
    // deactivate the buttons
    wv.curMode = -1;

    print_message("The server has terminated or network connection has been lost.","red");
}

function draw_colorbar(gl)
{
  if (wv.colorbar==undefined)
  {
    // remove the limits
    var textDiv_lo = document.getElementById("colorbar-text-lo");
    textDiv_lo.innerHTML = "";
    var textDiv_hi = document.getElementById("colorbar-text-hi");
    textDiv_hi.innerHTML = "";
    return;
  }

  gl.preserveDrawingBuffer = true;

  // Construct the identity as projection matrix and pass it in
  wv.mvpMatrix.load(wv.perspectiveMatrix);
  wv.mvpMatrix.setUniform(gl, wv.u_modelViewProjMatrixLoc, false);

  // get the colors sent from the server
  var ncol = wv.colorbar.length/3;

  var mv   = new J3DIMatrix4();
  var mVal = wv.mvMatrix.getAsArray();

  mVal[ 3] = 0.0;
  mVal[ 7] = 0.0;
  mVal[11] = 0.0;
  mv.load(mVal);
  mv.scale(1.0/wv.scale, 1.0/wv.scale, 1.0/wv.scale);
  mv.invert();
  mv.transpose();

  // define location of colorbar in space
  var x    = -1.3 * wv.width / wv.height;
  var y    =  0.05;
  var z    =  0.8;

  var w = 0.2;
  var h = 5*w;

  var ntri = 2*(ncol -1);
  var npts = 3*ntri;

  var dh = h/(ncol-1);

  var vertices = new Float32Array(3*npts);
  var colors = new Uint8Array(3*npts);
  var pt = 0;
  for (var i=0;i<ncol-1;i++)
  {
    // create the first triangle
    vertices[3*pt  ] = x;
    vertices[3*pt+1] = y +i*dh;
    vertices[3*pt+2] = z;
    for (var j=0;j<3;j++)
      colors[3*pt+j] = wv.colorbar[3*i+j];
    pt++;

    vertices[3*pt  ] = x +w;
    vertices[3*pt+1] = y +i*dh;
    vertices[3*pt+2] = z;
    for (var j=0;j<3;j++)
      colors[3*pt+j] = wv.colorbar[3*i+j];
    pt++;

    vertices[3*pt  ] = x +w;
    vertices[3*pt+1] = y +(i+1)*dh;
    vertices[3*pt+2] = z;
    for (var j=0;j<3;j++)
      colors[3*pt+j] = wv.colorbar[3*(i+1)+j];
    pt++;

    // create the second triangle
    vertices[3*pt  ] = x;
    vertices[3*pt+1] = y +i*dh;
    vertices[3*pt+2] = z;
    for (var j=0;j<3;j++)
      colors[3*pt+j] = wv.colorbar[3*i+j];
    pt++;

    vertices[3*pt  ] = x +w;
    vertices[3*pt+1] = y +(i+1)*dh;
    vertices[3*pt+2] = z;
    for (var j=0;j<3;j++)
      colors[3*pt+j] = wv.colorbar[3*(i+1)+j];
    pt++;

    vertices[3*pt  ] = x;
    vertices[3*pt+1] = y +(i+1)*dh;
    vertices[3*pt+2] = z;
    for (var j=0;j<3;j++)
      colors[3*pt+j] = wv.colorbar[3*(i+1)+j];
    pt++;
  }

  // draw the triangles defining the colorbar
  gl.disableVertexAttribArray(2);
  gl.disableVertexAttribArray(1);

  gl.uniform1f(wv.u_wLightLoc, 0.0);

  // create the vertex buffer
  var buffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
  gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);
  gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
  gl.enableVertexAttribArray(0);

  //gl.uniform1f(wv.u_wNormalLoc, 0.0);

  // create the color buffer
  var cbuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, cbuf);
  gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STATIC_DRAW);
  gl.vertexAttribPointer(1, 3, gl.UNSIGNED_BYTE, false, 0, 0);
  gl.enableVertexAttribArray(1);

  gl.drawArrays(gl.TRIANGLES, 0, npts );
  gl.deleteBuffer(buffer);
  gl.deleteBuffer(cbuf);
  gl.uniform1f(wv.u_wLightLoc, 1.0);

  wv.checkGLError(gl, "draw_colorbar");

  var fs = 16;
  var xlim_lo = [x+1.1*w,y+0.05*h,z];
  var xlim_hi = [x+1.1*w,y+h,z];

  var xlims = [xlim_lo,xlim_hi];
  var divs = ["colorbar-text-lo","colorbar-text-hi"];

  var M = wv.mvpMatrix.getAsArray();

  for (var i=0;i<xlims.length;i++)
  {
    var X = new J3DIVector3(xlims[i][0],xlims[i][1],xlims[i][2]);
    var b = new J3DIVector3();

    // transpose(mvp)*(xyz)
    b[0] = M[0]*X[0] +M[4]*X[1] +M[8]* X[2] +M[12]*1;
    b[1] = M[1]*X[0] +M[5]*X[1] +M[9]* X[2] +M[13]*1;
    b[2] = M[2]*X[0] +M[6]*X[1] +M[10]*X[2] +M[14]*1;
    b[3] = M[3]*X[0] +M[7]*X[1] +M[11]*X[2] +M[15]*1;

    // normalize like the gpu
    b[0] /= b[3];
    b[1] /= b[3];
    b[2] /= b[3];

    // compute the pixel coordinates
    var px = Math.round((( b[0] + 1 ) / 2.0) * wv.width );
    var py = Math.round((( 1 -b[1]) / 2.0) * wv.height );

    var textDiv = document.getElementById(divs[i]);
    textDiv.style.left = px +"px";
    textDiv.style.top  = py +"px";
    textDiv.style.fontSize = fs.toString()+"px";
    if (document.getElementById("background-checkbox").checked)
      textDiv.style.color = "white";
    else
      textDiv.style.color = "black";
    textDiv.innerHTML = wv.colorbarlims[i].toString();
  }
}

//
// callback used to put axes on the canvas (called by webplot.html)
//
function wvUpdateCanvas(gl)
{
    gl.preserveDrawingBuffer = true;

    // Construct the identity as projection matrix and pass it in
    wv.mvpMatrix.load(wv.perspectiveMatrix);
    wv.mvpMatrix.setUniform(gl, wv.u_modelViewProjMatrixLoc, false);

    var mv   = new J3DIMatrix4();
    var mVal = wv.mvMatrix.getAsArray();

    mVal[ 3] = 0.0;
    mVal[ 7] = 0.0;
    mVal[11] = 0.0;
    mv.load(mVal);
    mv.scale(1.0/wv.scale, 1.0/wv.scale, 1.0/wv.scale);
    mv.invert();
    mv.transpose();

    // define location of axes in space
    var x    = -1.25 * wv.width / wv.height;
    var y    = -1.25;
    var z    =  0.8;
    var alen =  0.25;     // length of axes

    // set up coordinates for axes
    mVal = mv.getAsArray();

    var vertices = new Float32Array(66);
    vertices[ 0] = x;
    vertices[ 1] = y;
    vertices[ 2] = z;
    vertices[ 3] = x + alen*(    mVal[ 0]             );
    vertices[ 4] = y + alen*(    mVal[ 1]             );
    vertices[ 5] = z + alen*(    mVal[ 2]             );
    var drawX = !wv.mode4d;
    if (drawX)
    {
      vertices[ 6] = x + alen*(1.1*mVal[ 0]+0.1*mVal[ 4]);
      vertices[ 7] = y + alen*(1.1*mVal[ 1]+0.1*mVal[ 5]);
      vertices[ 8] = z + alen*(1.1*mVal[ 2]+0.1*mVal[ 6]);
      vertices[ 9] = x + alen*(1.3*mVal[ 0]-0.1*mVal[ 4]);
      vertices[10] = y + alen*(1.3*mVal[ 1]-0.1*mVal[ 5]);
      vertices[11] = z + alen*(1.3*mVal[ 2]-0.1*mVal[ 6]);
      vertices[12] = x + alen*(1.1*mVal[ 0]-0.1*mVal[ 4]);
      vertices[13] = y + alen*(1.1*mVal[ 1]-0.1*mVal[ 5]);
      vertices[14] = z + alen*(1.1*mVal[ 2]-0.1*mVal[ 6]);
      vertices[15] = x + alen*(1.3*mVal[ 0]+0.1*mVal[ 4]);
      vertices[16] = y + alen*(1.3*mVal[ 1]+0.1*mVal[ 5]);
      vertices[17] = z + alen*(1.3*mVal[ 2]+0.1*mVal[ 6]);
    }
    else
    {
      // draw a U
    }

    vertices[18] = x;
    vertices[19] = y;
    vertices[20] = z;
    vertices[21] = x + alen*(    mVal[ 4]             );
    vertices[22] = y + alen*(    mVal[ 5]             );
    vertices[23] = z + alen*(    mVal[ 6]             );
    if (drawX)
    {
      vertices[24] = x + alen*(1.1*mVal[ 4]+0.1*mVal[ 8]);
      vertices[25] = y + alen*(1.1*mVal[ 5]+0.1*mVal[ 9]);
      vertices[26] = z + alen*(1.1*mVal[ 6]+0.1*mVal[10]);
      vertices[27] = x + alen*(1.2*mVal[ 4]             );
      vertices[28] = y + alen*(1.2*mVal[ 5]             );
      vertices[29] = z + alen*(1.2*mVal[ 6]             );
      vertices[30] = x + alen*(1.3*mVal[ 4]+0.1*mVal[ 8]);
      vertices[31] = y + alen*(1.3*mVal[ 5]+0.1*mVal[ 9]);
      vertices[32] = z + alen*(1.3*mVal[ 6]+0.1*mVal[10]);
      vertices[33] = x + alen*(1.2*mVal[ 4]             );
      vertices[34] = y + alen*(1.2*mVal[ 5]             );
      vertices[35] = z + alen*(1.2*mVal[ 6]             );
      vertices[36] = x + alen*(1.2*mVal[ 4]             );
      vertices[37] = y + alen*(1.2*mVal[ 5]             );
      vertices[38] = z + alen*(1.2*mVal[ 6]             );
      vertices[39] = x + alen*(1.2*mVal[ 4]-0.1*mVal[ 8]);
      vertices[40] = y + alen*(1.2*mVal[ 5]-0.1*mVal[ 9]);
      vertices[41] = z + alen*(1.2*mVal[ 6]-0.1*mVal[10]);
    }
    else
    {
      // draw a V
    }

    vertices[42] = x;
    vertices[43] = y;
    vertices[44] = z;
    vertices[45] = x + alen*(    mVal[ 8]             );
    vertices[46] = y + alen*(    mVal[ 9]             );
    vertices[47] = z + alen*(    mVal[10]             );
    if (drawX)
    {
      vertices[48] = x + alen*(1.1*mVal[ 8]+0.1*mVal[ 0]);
      vertices[49] = y + alen*(1.1*mVal[ 9]+0.1*mVal[ 1]);
      vertices[50] = z + alen*(1.1*mVal[10]+0.1*mVal[ 2]);
      vertices[51] = x + alen*(1.3*mVal[ 8]+0.1*mVal[ 0]);
      vertices[52] = y + alen*(1.3*mVal[ 9]+0.1*mVal[ 1]);
      vertices[53] = z + alen*(1.3*mVal[10]+0.1*mVal[ 2]);
      vertices[54] = x + alen*(1.3*mVal[ 8]+0.1*mVal[ 0]);
      vertices[55] = y + alen*(1.3*mVal[ 9]+0.1*mVal[ 1]);
      vertices[56] = z + alen*(1.3*mVal[10]+0.1*mVal[ 2]);
      vertices[57] = x + alen*(1.1*mVal[ 8]-0.1*mVal[ 0]);
      vertices[58] = y + alen*(1.1*mVal[ 9]-0.1*mVal[ 1]);
      vertices[59] = z + alen*(1.1*mVal[10]-0.1*mVal[ 2]);
      vertices[60] = x + alen*(1.1*mVal[ 8]-0.1*mVal[ 0]);
      vertices[61] = y + alen*(1.1*mVal[ 9]-0.1*mVal[ 1]);
      vertices[62] = z + alen*(1.1*mVal[10]-0.1*mVal[ 2]);
      vertices[63] = x + alen*(1.3*mVal[ 8]-0.1*mVal[ 0]);
      vertices[64] = y + alen*(1.3*mVal[ 9]-0.1*mVal[ 1]);
      vertices[65] = z + alen*(1.3*mVal[10]-0.1*mVal[ 2]);
    }
    else
    {
      // draw a W
    }

    // set up colors for the axes
    var colors = new Uint8Array(66);
    colors[ 0] = 255;   colors[ 1] =   0;   colors[ 2] =   0;
    colors[ 3] = 255;   colors[ 4] =   0;   colors[ 5] =   0;
    colors[ 6] = 255;   colors[ 7] =   0;   colors[ 8] =   0;
    colors[ 9] = 255;   colors[10] =   0;   colors[11] =   0;
    colors[12] = 255;   colors[13] =   0;   colors[14] =   0;
    colors[15] = 255;   colors[16] =   0;   colors[17] =   0;
    colors[18] =   0;   colors[19] = 255;   colors[20] =   0;
    colors[21] =   0;   colors[22] = 255;   colors[23] =   0;
    colors[24] =   0;   colors[25] = 255;   colors[26] =   0;
    colors[27] =   0;   colors[28] = 255;   colors[29] =   0;
    colors[30] =   0;   colors[31] = 255;   colors[32] =   0;
    colors[33] =   0;   colors[34] = 255;   colors[35] =   0;
    colors[36] =   0;   colors[37] = 255;   colors[38] =   0;
    colors[39] =   0;   colors[40] = 255;   colors[41] =   0;
    colors[42] =   0;   colors[43] =   0;   colors[44] = 255;
    colors[45] =   0;   colors[46] =   0;   colors[47] = 255;
    colors[48] =   0;   colors[49] =   0;   colors[50] = 255;
    colors[51] =   0;   colors[52] =   0;   colors[53] = 255;
    colors[54] =   0;   colors[55] =   0;   colors[56] = 255;
    colors[57] =   0;   colors[58] =   0;   colors[59] = 255;
    colors[60] =   0;   colors[61] =   0;   colors[62] = 255;
    colors[63] =   0;   colors[64] =   0;   colors[65] = 255;

    // draw the axes
    if (gl.lineWidth) {
        gl.lineWidth(3);
    }
    gl.disableVertexAttribArray(2);
    gl.uniform1f(wv.u_wLightLoc, 0.0);

    var buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(0);

    var cbuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, cbuf);
    gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STATIC_DRAW);
    gl.vertexAttribPointer(1, 3, gl.UNSIGNED_BYTE, false, 0, 0);
    gl.enableVertexAttribArray(1);

    gl.drawArrays(gl.LINES, 0, 22);
    gl.deleteBuffer(buffer);
    gl.deleteBuffer(cbuf);
    gl.uniform1f(wv.u_wLightLoc, 1.0);

    draw_colorbar(gl);

    // request that the scene be updated after drawing the colorbar
    wv.sceneUpd = 1;
}

function sleep(milliseconds) {
  var start = new Date().getTime();
  for (var i = 0; i < 1e7; i++) {
    if ((new Date().getTime() - start) > milliseconds){
      break;
    }
  }
}

//
// callback when "homeButton" is pressed (called by webplot.html)
//
function cmdHome()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    }
}


//
// callback when "leftButton" is pressed (called by webplot.html)
//
function cmdLeft()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.mvMatrix.rotate(+90, 0,1,0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    }
}


//
// callback when "riteButton" is pressed (called by webplot.html)
//
function cmdRite()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.mvMatrix.rotate(-90, 0,1,0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    }
}


//
// callback when "botmButton" is pressed (called by webplot.html)
//
function cmdBotm()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.mvMatrix.rotate(-90, 1,0,0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    }
}


//
// callback when "topButton" is pressed (called by webplot.html)
//
function cmdTop()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.makeIdentity();
        wv.mvMatrix.rotate(+90, 1,0,0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale    = 1;
        wv.sceneUpd = 1;
    }
}


//
// callback when "inButton" is pressed (called by webplot.html)
//
function cmdIn()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.load(wv.uiMatrix);
        wv.mvMatrix.scale(2.0, 2.0, 2.0);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale   *= 2.0;
        wv.sceneUpd = 1;
    }
}


//
// callback when "outButton" is pressed (called by webplot.html)
//
function cmdOut()
{
    if (wv.curMode == 0) {
        wv.mvMatrix.load(wv.uiMatrix);
        wv.mvMatrix.scale(0.5, 0.5, 0.5);
        wv.uiMatrix.load(wv.mvMatrix);
        wv.mvMatrix.makeIdentity();
        wv.scale   *= 0.5;
        wv.sceneUpd = 1;
    }
}


//
// callback when any mouse is pressed in canvas (when wv.curMode==0)
//
function getMouseDown(e)
{
    if (!e) var e = event;

    wv.startX   =  e.clientX - wv.offLeft             - 1;
    wv.startY   = -e.clientY + wv.offTop  + wv.height + 1;

    wv.dragging = true;
    wv.button   = e.button;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
}


//
// callback when the mouse moves in canvas (when wv.curMode==0)
//
function getMouseMove(e)
{
    if (!e) var e = event;

    wv.cursorX  =  e.clientX - wv.offLeft             - 1;
    wv.cursorY  = -e.clientY + wv.offTop  + wv.height + 1;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
}


//
// callback when the mouse is released in canvas (when wv.curMode==0)
//
function getMouseUp(e)
{
    wv.dragging = false;
}


//
// callback when the mouse wheel is rolled in canvas (when wv.curMode==0)
//
function getMouseRoll(e)
{
    if (e) {

        // zoom in
        if        (e.deltaY > 0) {
            wv.mvMatrix.scale(1.1, 1.1, 1.1);
            wv.scale *= 1.1;
            wv.sceneUpd = 1;

        // zoom out
        } else if (e.deltaY < 0) {
            wv.mvMatrix.scale(0.9, 0.9, 0.9);
            wv.scale *= 0.9;
            wv.sceneUpd = 1;
        }
    }
}


//
// callback when the mouse leaves the canvas (when wv.curMode==0)
//
function mouseLeftCanvas(e)
{
    if (wv.dragging) {
        wv.dragging = false;
    }
}

//
// callback when a key is pressed
//
function getKeyPress(e)
{
    if (wv.curMode == 0) {
        wv.keyPress = e.charCode;
        wv.keyCode  = e.keyCode;

                        wv.modifier  = 0;
        if (e.shiftKey) wv.modifier |= 1;
        if (e.altKey  ) wv.modifier |= 2;
        if (e.ctrlKey ) wv.modifier |= 4;

    }

    return true;
}


//
// callback when an arrow... or shift is pressed (needed for Chrome)
//
function getKeyDown(e)
{
    if (e.charCode == 0) {
        if (e.keyCode == 33 ||          // PgUp
            e.keyCode == 34 ||          // PgDn
            e.keyCode == 36 ||          // Home
            e.keyCode == 37 ||          // Left
            e.keyCode == 38 ||          // Up
            e.keyCode == 39 ||          // Right
            e.keyCode == 40   ) {       // Down
            wv.keyCode  = e.keyCode;
            wv.keyPress = 0;
        } else if (e.keyCode == 16) {   // Shift
            wv.modifier = 1;
        }
    }

}

//
// callback when a shift is released (needed for Chrome)
//
function getKeyUp(e)
{

    if (e.charCode == 0 && e.keyCode == 16) {
        wv.modifier = 0;
    }
}

//
// print an object and its contents
//
function printObject(obj)
{
    var out = '';

    for (var p in obj) {
        out += p + ': ' + obj[p] + '\n';
    }

    alert(out);
}


//
// simple text formatter patterned after C's sprintf
//
function sprintf()
{

    // if no arguments, return an empty string
    var narg = arguments.length;
    if (narg == 0) {
        return "";
    }

    // otherwise, build output from input
    var format = arguments[0];
    var answer = "";
    var iarg   = 1;

    while (1) {
        var indx = format.indexOf("%");
        if (indx >= 0) {
            answer += format.substring(0, indx);
            if (iarg < narg) {
                answer += arguments[iarg++];
            } else {
                answer += "%";
            }
            format  = format.substring(indx+1);
        } else {
            answer += format;
            break;
        }
    }

    return answer;
}
