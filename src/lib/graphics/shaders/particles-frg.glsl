//#version 410

layout( location = 0 ) out vec4 fragColor;

in float g_Color;

void main() {
  if (g_Color < 2.0) {
    //discard;
    fragColor = vec4(0,0,1,0.5);
  }
  else
    fragColor = vec4(1,0,0,0.5);
}
