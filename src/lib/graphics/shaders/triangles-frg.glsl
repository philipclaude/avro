#version 410

layout( location = 0 ) out vec4 fragColor;

in vec2 v_Parameter;
in vec3 g_Normal;
in vec3 g_Position;
in vec3 x_Position;

uniform samplerBuffer solution;
uniform samplerBuffer colormap;

uniform vec3 constant_color;
uniform int use_constant_color;
uniform int have_tessellation_shader;
uniform int u_lighting;
uniform float u_alpha;
uniform int u_clip;
uniform vec3 u_clip_center;
uniform vec3 u_clip_normal;

flat in float g_clip;

// TODO: make these uniforms
const int ncolor = 256;

uniform float u_umin;
uniform float u_umax;

void
get_color( float u , out vec3 color ) {

  float umin = u_umin;
  float umax = u_umax;

    int indx = int(ncolor*(u - umin)/(umax - umin));

    if (indx < 0) indx = 0;
    if (indx > 255) indx = 255;

    float r0 = texelFetch( colormap , 3*(indx) + 0 ).x;
    float g0 = texelFetch( colormap , 3*(indx) + 1 ).x;
    float b0 = texelFetch( colormap , 3*(indx) + 2 ).x;

    color = vec3(r0,g0,b0);
}

void
shading( in vec3 l , in vec3 n , in vec3 color , out vec3 color_out ) {

  float diffuse = max(0.0,dot(l,n));
  float phong = 128.0;
  float specular = pow(max(0.0,dot(-reflect(l,n),n)),phong);

  vec3 cd = color * diffuse;
  vec3 cs = vec3(0.2) * specular;
  vec3 ca = vec3(0.2);
  color_out = ca + cd + cs;
}

void main() {

  if (u_clip == 2 && g_clip < 0.0) discard;
  if (u_clip > 0) {

    float p = dot(x_Position - u_clip_center,u_clip_normal);
    if (p < 0.0)
      discard;
  }

  if (use_constant_color > 0) {
    vec3 color_out;
    if (have_tessellation_shader > 0 && u_lighting > 0)
      shading( -normalize(g_Position) , normalize(g_Normal) , constant_color , color_out );
    else
      color_out = constant_color; // no shading because there are no normals
    fragColor = vec4(color_out,u_alpha);
    return;
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

  float phi0 =  s*(-1.1E1/2.0)-t*(1.1E1/2.0)+s*t*1.8E1-s*(t*t)*(2.7E1/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*9.0-
              (s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*t*t)*(9.0/2.0)+1.0;
  float phi1 =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);
  float phi2 =  t-(t*t)*(9.0/2.0)+(t*t*t)*(9.0/2.0);
  float phi3 =  s*t*(-9.0/2.0)+(s*s)*t*(2.7E1/2.0);
  float phi4 =  s*t*(-9.0/2.0)+s*(t*t)*(2.7E1/2.0);
  float phi5 =  t*(-9.0/2.0)+s*t*(9.0/2.0)-s*(t*t)*(2.7E1/2.0)+(t*t)*1.8E1-(t*t*t)*(2.7E1/2.0);
  float phi6 =  t*9.0-s*t*(4.5E1/2.0)+s*(t*t)*2.7E1+(s*s)*t*(2.7E1/2.0)-(t*t)*(4.5E1/2.0)+(t*t*t)*(2.7E1/2.0);
  float phi7 =  s*9.0-s*t*(4.5E1/2.0)+s*(t*t)*(2.7E1/2.0)+(s*s)*t*2.7E1-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0);
  float phi8 =  s*(-9.0/2.0)+s*t*(9.0/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*1.8E1-(s*s*s)*(2.7E1/2.0);
  float phi9 =  s*t*2.7E1-s*(t*t)*2.7E1-(s*s)*t*2.7E1;

  f = f0*phi0 + f1*phi1 + f2*phi2 + f3*phi3 + f4*phi4 + f5*phi5 + f6*phi6 + f7*phi7 + f8*phi8 + f9*phi9;
  get_color(f,color);

  #else
  color = vec3(0.8,0.4,0.2);
  #endif

  vec3 color_out = color;
  if (u_lighting > 0)
    shading( -normalize(g_Position) , normalize(g_Normal) , color , color_out );

  fragColor = vec4(color_out,u_alpha);
}
