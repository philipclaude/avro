#version 410

uniform int u_clip;
uniform vec3 u_clip_center;
uniform vec3 u_clip_normal;

in vec3 x_Position;
flat in float g_clip;

layout( location = 0 ) out vec4 fragColor;

void main() {

  if (u_clip == 2 && g_clip < 0.0) discard;
  if (u_clip > 0) {
    float p = dot(x_Position - u_clip_center,u_clip_normal);
    if (p < 0.0)
      discard;
  }

  // draw edges in black
  fragColor = vec4(0,0,0,1);
}
