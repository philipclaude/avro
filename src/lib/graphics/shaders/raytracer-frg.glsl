//#version 410

layout( location = 0 ) out vec4 fragColor;

in vec2 v_Parameter;

uniform samplerBuffer pixels;

uniform int u_width;
uniform int u_height;

void main() {

  int j = int( gl_FragCoord.x - 0.5)/2;
  int i = int( gl_FragCoord.y - 0.5)/2;

  j = int(v_Parameter.x * u_width);
  i = int(v_Parameter.y * u_height);

  int idx = i * u_width + j;
  vec3 color = texelFetch( pixels , idx ).xyz;

  fragColor = vec4(color,1.0);
}
