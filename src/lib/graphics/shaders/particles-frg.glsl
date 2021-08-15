//#version 410

layout( location = 0 ) out vec4 fragColor;

in float v_Color;

void main() {
  if (v_Color < 2.0)
    fragColor = vec4(0,0,1,0.5);
  else
    fragColor = vec4(1,0,0,0.5);
}
