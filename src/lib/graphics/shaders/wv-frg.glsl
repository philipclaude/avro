#version 330

precision mediump float;
precision mediump int;

// uniforms
uniform float bColor;                       // backface color switch
uniform int   picking;                      // picking flag
uniform int   vbonum;                       // vbo number

// inputs from the vertex shader
in float z_Screen;
in vec4  v_Color;
in vec4  v_bColor;

// outputs
out vec4 FragColor;

void main()
{
   if (picking == 0)
   {
     FragColor = v_Color;
     if ((bColor != 0.0) && (gl_FrontFacing==false))
       FragColor = v_bColor;
   }
   else if (picking < 0)
   {
     FragColor = vec4(0.0, 0.0, 0.0, 0.0);
     // save away x/100, y/100, z or 1/w value
     float              compnt = gl_FragCoord.x/100.;    // pixel
     if (picking == -2) compnt = gl_FragCoord.y/100.;    // pixel
     if (picking == -3) compnt = gl_FragCoord.w*z_Screen;
     if (picking == -4) compnt = gl_FragCoord.w;
     int pmsign = 0;
     if (compnt < 0.0)
     {
       compnt = -compnt;
       pmsign = 1;
     }
     if (compnt > 127.0) return;
     FragColor.r = floor(compnt);
     if (pmsign == 1) FragColor.r += 128.0;
     compnt         = fract(compnt)*256.0;
     FragColor.g = floor(compnt);
     compnt         = fract(compnt)*256.0;
     FragColor.b = floor(compnt);
     FragColor.a = floor(fract(compnt)*256.0);
     if (FragColor.a == 256.0)
     {
       FragColor.a  = 0.0;
       FragColor.b += 1.0;
     }
     if (FragColor.b >= 256.0)
     {
       FragColor.b -= 256.0;
       FragColor.g += 1.0;
     }
     if (FragColor.g >= 256.0)
     {
       FragColor.g -= 256.0;
       FragColor.r += 1.0;
     }
     if (FragColor.a == 0.0) FragColor.a = 1.0;
     FragColor /= 255.;
   }
   else
   {
     int high       = vbonum/256;
     FragColor.r = float(high)/255.;
     FragColor.g = float(vbonum - high*256)/255.;
     FragColor.b = 0.0;
     FragColor.a = 0.0;
//        high           = gl_PrimitiveID/256;
//        FragColor.b = float(high)/255.;
//        FragColor.a = float(gl_PrimitiveID - high*256)/255.;
   }
}
