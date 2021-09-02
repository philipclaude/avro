//#version 410

layout( location = 0 ) out vec4 fragColor;

in vec2 v_Parameter;

uniform samplerBuffer pixels;

uniform int u_width;
uniform int u_height;

void main() {

  float j = v_Parameter.x * u_width;
  float i = v_Parameter.y * u_height;

  int idx = int(i * u_width + j);
  vec3 color = texelFetch( pixels , idx ).xyz;

  fragColor = vec4(v_Parameter.x,0,0,1.0);
  fragColor = vec4(color,1.0);
}
