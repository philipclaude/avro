#version 410


layout( location = 0 ) out vec4 fragColor;

in vec2 v_Parameter;
in vec3 g_Normal;
in vec3 g_Position;

uniform samplerBuffer solution;
uniform samplerBuffer colormap;

uniform vec3 constant_color;
uniform int use_constant_color;

// TODO: make these uniforms
const int ncolor = 256;
const float umin = -1;
const float umax =  1;

void
get_color( float u , out vec3 color ) {

    int indx = int(ncolor*(u - umin)/(umax - umin));

    float r0 = texelFetch( colormap , 3*(indx) + 0 ).x;
    float g0 = texelFetch( colormap , 3*(indx) + 1 ).x;
    float b0 = texelFetch( colormap , 3*(indx) + 2 ).x;

    color = vec3(r0,g0,b0);
}

void main() {

  if (use_constant_color > 0) {
    //fragColor = vec4(constant_color,1.0);
    //fragColor = vec4( abs(g_Normal),1.0);
    //return;
  }

  float s = v_Parameter.x;
  float t = v_Parameter.y;

  float f = 0.0;
  vec3 color;

  #if SOLUTION_ORDER == 0

  int idx = gl_PrimitiveID*1;
  f = texelFetch( solution , idx + 0 ).x;
  get_color(f,color);

  #elif SOLUTION_ORDER == 1

  int idx  = gl_PrimitiveID*3;

  float f0 = texelFetch( solution , idx + 0 ).x;
  float f1 = texelFetch( solution , idx + 1 ).x;
  float f2 = texelFetch( solution , idx + 2 ).x;

  f = (1 - s - t)*f0 + s*f1 + t*f2;
  get_color(f,color);

  #elif SOLUTION_ORDER == 2

  int idx  = gl_PrimitiveID*6;

  float f0 = texelFetch( solution , idx + 0 ).x;
  float f1 = texelFetch( solution , idx + 1 ).x;
  float f2 = texelFetch( solution , idx + 2 ).x;
  float f3 = texelFetch( solution , idx + 3 ).x;
  float f4 = texelFetch( solution , idx + 4 ).x;
  float f5 = texelFetch( solution , idx + 5 ).x;

  float phi0 =  s*-3.0-t*3.0+s*t*4.0+(s*s)*2.0+(t*t)*2.0+1.0;
  float phi1 =  -s+(s*s)*2.0;
  float phi2 =  -t+(t*t)*2.0;
  float phi3 =  s*t*4.0;
  float phi4 =  t*(s+t-1.0)*-4.0;
  float phi5 =  -s*(s*4.0+t*4.0-4.0);
  f = f0*phi0 + f1*phi1 + f2*phi2 + f3*phi3 + f4*phi4 + f5*phi5;
  get_color(f,color);

  #elif SOLUTION_ORDER == 3

  int idx  = gl_PrimitiveID*10;

  float f0 = texelFetch( solution , idx + 0 ).x;
  float f1 = texelFetch( solution , idx + 1 ).x;
  float f2 = texelFetch( solution , idx + 2 ).x;
  float f3 = texelFetch( solution , idx + 3 ).x;
  float f4 = texelFetch( solution , idx + 4 ).x;
  float f5 = texelFetch( solution , idx + 5 ).x;
  float f6 = texelFetch( solution , idx + 6 ).x;
  float f7 = texelFetch( solution , idx + 7 ).x;
  float f8 = texelFetch( solution , idx + 8 ).x;
  float f9 = texelFetch( solution , idx + 9 ).x;

  float u[10];
  u[0] = f0; u[1] = f1; u[2] = f2; u[3] = f3; u[4] = f4; u[5] = f5; u[6] = f6; u[7] = f7; u[8] = f8; u[9] = f9;
  float phi[10];
  phi[0] =  s*(-1.1E1/2.0)-t*(1.1E1/2.0)+s*t*1.8E1-s*(t*t)*(2.7E1/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*9.0-
              (s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*t*t)*(9.0/2.0)+1.0;

  phi[1] =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);

  phi[2] =  t-(t*t)*(9.0/2.0)+(t*t*t)*(9.0/2.0);

  phi[3] =  s*t*(-9.0/2.0)+(s*s)*t*(2.7E1/2.0);

  phi[4] =  s*t*(-9.0/2.0)+s*(t*t)*(2.7E1/2.0);

  phi[5] =  t*(-9.0/2.0)+s*t*(9.0/2.0)-s*(t*t)*(2.7E1/2.0)+(t*t)*1.8E1-(t*t*t)*(2.7E1/2.0);

  phi[6] =  t*9.0-s*t*(4.5E1/2.0)+s*(t*t)*2.7E1+(s*s)*t*(2.7E1/2.0)-(t*t)*(4.5E1/2.0)+(t*t*t)*(2.7E1/2.0);

  phi[7] =  s*9.0-s*t*(4.5E1/2.0)+s*(t*t)*(2.7E1/2.0)+(s*s)*t*2.7E1-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0);

  phi[8] =  s*(-9.0/2.0)+s*t*(9.0/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*1.8E1-(s*s*s)*(2.7E1/2.0);

  phi[9] =  s*t*2.7E1-s*(t*t)*2.7E1-(s*s)*t*2.7E1;

  f = u[0]*phi[0] + u[1]*phi[1] + u[2]*phi[2] + u[3]*phi[3] + u[4]*phi[4] + u[5]*phi[5] + u[6]*phi[6] + u[7]*phi[7] + u[8]*phi[8] + u[9]*phi[9];
  get_color(f,color);

  #else
  color = vec3(0.8,0.8,0.2);
  #endif

  #if 1
  // vector from surface point to camera (where light is)
  vec3 direction = -normalize(g_Position);
  vec3 normal = normalize(g_Normal);

  float diffuse = max(0.0,dot(direction,normal));
  float phong = 128.0;
  float specular = pow(max(0.0,dot(-reflect(direction,normal),normal)),phong);

  //color = constant_color;
  //color = vec3(0.4);
  vec3 cd = color * diffuse;
  vec3 cs = vec3(0.2) * specular;
  vec3 ca = vec3(0.2);
  color = ca + cd + cs;
  #endif

  fragColor = vec4(color,1.0);
}
