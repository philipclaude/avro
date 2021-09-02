//#version 330 core

layout(points) in;
layout(triangle_strip, max_vertices = 6) out;

out vec2 v_Parameter;

void main()
{
  // first triangle
  gl_Position = vec4( -1.0 , -1.0 , 0.0 , 1.0 );
  v_Parameter = vec2( 0.0 , 0.0 );
  EmitVertex();

  gl_Position = vec4( 1.0 , -1.0 , 0.0 , 1.0 );
  v_Parameter = vec2( 1.0 , 0.0 );
  EmitVertex();

  gl_Position = vec4( 1.0 , 1.0 , 0.0 , 1.0 );
  v_Parameter = vec2( 1.0 , 1.0 );
  EmitVertex();

  EndPrimitive();

  // second triangle
  gl_Position = vec4( -1.0 , -1.0 , 0.0 , 1.0 );
  v_Parameter = vec2( 0.0 , 0.0 );
  EmitVertex();

  gl_Position = vec4( 1.0 , 1.0 , 0.0 , 1.0 );
  v_Parameter = vec2( 1.0 , 1.0 );
  EmitVertex();

  gl_Position = vec4( -1.0 , 1.0 , 0.0 , 1.0 );
  v_Parameter = vec2( 0.0 , 1.0 );
  EmitVertex();

  EndPrimitive();

}
