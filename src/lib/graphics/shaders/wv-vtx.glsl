#version 330

precision mediump float;
precision mediump int;

// uniforms
uniform mat4   u_modelViewMatrix;           // not currently used
uniform mat4   MVP;                         // model-view-proj matrix
uniform mat4   u_normalMatrix;
uniform vec3   lightDir;
uniform vec3   conNormal;                   // constant normal
uniform vec3   conColor;                    // constant color
uniform vec3   bacColor;                    // back face color
uniform float  wAmbient;                    // Ambient light weight
uniform float  wColor;                      // Constant color switch
uniform float  bColor;                      // Backface color switch
uniform float  wNormal;                     // Constant normal switch
uniform float  wLight;                      // lighting switch
uniform float  xpar;                        // transparency factor
uniform float  pointSize;                   // point size in pixels
uniform int    picking;                     // picking flag

// inputs
layout (location = 0 ) in vec4 vPosition;
layout (location = 1 ) in vec4 vColor;
layout (location = 2 ) in vec3 vNormal;

// outputs
out float z_Screen;
out vec4 v_Color;
out vec4 v_bColor;

void main()
{
   // set the pixel position
   gl_Position = MVP*vPosition;
   z_Screen    = gl_Position[2];

   // set the point size
   if (wLight <= 0.0) gl_PointSize = pointSize;

   // return if picking is off
   if (picking != 0) return;

   // assumes that colors are coming in as unsigned bytes
   vec4 color = vColor/255.0;

   if (wLight <= 0.0)
   {
     // no lighting
     v_Color  = color*wColor + vec4(conColor,1)*(1.0-wColor);
     v_bColor = v_Color;
   }
   else
   {
     // setup bi-directional lighting
     //   a simple ambient/diffuse lighting model is used with:
     //      - single 'white' source & no 'material' color
     //      - linear mixture of ambient & diffuse based on weight
     vec3 lDirection = normalize(lightDir);
     vec3 norm       = vNormal*wNormal + conNormal*(1.0-wNormal);
     vec3 normal     = normalize(u_normalMatrix * vec4(norm,1)).xyz;
     float dp        = abs(dot(normal, lDirection));

     // make the color to be rendered
     color           = color*wColor + vec4(conColor,1)*(1.0-wColor);
     v_Color         = color*dp + color*wAmbient;
     v_bColor        = v_Color;

     // are we coloring the backface?
     if (bColor != 0.0)
     {
         color       = vec4(bacColor,1);
         v_bColor    = color*dp + color*wAmbient;
     }
   }
   v_Color.a  = xpar;
   v_bColor.a = xpar;
}
