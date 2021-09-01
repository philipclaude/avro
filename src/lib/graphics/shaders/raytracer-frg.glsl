//#version 410

layout( location = 0 ) out vec4 fragColor;

in vec2 v_Parameter;

uniform samplerBuffer pixels;

uniform int u_width;
uniform int u_height;

void main() {

  float u = v_Parameter.x * u_width;
  float v = v_Parameter.y * u_height;

  int idx = int(v * u_width + u);
  vec3 color = texelFetch( pixels , idx ).xyz;

  //fragColor = vec4(u,v,0.,1.);
  fragColor = vec4(color,1.0);
}
