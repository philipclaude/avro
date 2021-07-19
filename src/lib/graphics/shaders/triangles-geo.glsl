//#version 410 core

layout (triangles) in;
uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_NormalMatrix;
uniform mat4 u_ModelViewMatrix;
uniform int u_clip;
uniform vec3 u_clip_center;
uniform vec3 u_clip_normal;

in vec3 v_Position[]; // receive the vertex coordinates from the vertex shader or tessellation shader
in vec3 v_Normal[];

out vec2 v_Parameter; // pass the parameter coordinates to the fragment shader
out vec3 g_Position;
out vec3 g_Normal;
out vec3 x_Position;

in vec2 v_ParameterTess[];

flat out float g_clip;
flat in float v_clip[];

layout (triangle_strip , max_vertices = 3) out;

float
visible( in vec3 x ) {
  float result = 1.0;
  if (dot( x - u_clip_center , u_clip_normal ) < 0.0)
    result = -1.0;
  return result;
}

void main() {

  //vec3 normal = mat3(u_NormalMatrix) * cross( v_Position[1] - v_Position[0] , v_Position[2] - v_Position[0] );

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[0],1.0);
  v_Parameter = v_ParameterTess[0];
  g_Position  = (u_ModelViewMatrix * vec4(v_Position[0],1.0)).xyz;
  g_Normal    = mat3(u_NormalMatrix) * v_Normal[0];
  x_Position  = v_Position[0];
  g_clip      = (u_clip == 2) ? v_clip[0] : (u_clip == 1) ? visible( v_Position[0] ) : 0;
  #if WITH_TESSELLATION == 0
  v_Parameter = vec2(0,0);
  #endif
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[1],1.0);
  v_Parameter = v_ParameterTess[1];
  g_Position  = (u_ModelViewMatrix * vec4(v_Position[1],1.0)).xyz;
  g_Normal    = mat3(u_NormalMatrix) * v_Normal[1];
  x_Position  = v_Position[1];
  g_clip      = (u_clip == 2) ? v_clip[1] : (u_clip == 1) ? visible( v_Position[1] ) : 0;
  #if WITH_TESSELLATION == 0
  v_Parameter = vec2(1,0);
  #endif
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[2],1.0);
  v_Parameter = v_ParameterTess[2];
  g_Position  = (u_ModelViewMatrix * vec4(v_Position[2],1.0)).xyz;
  g_Normal    = mat3(u_NormalMatrix) * v_Normal[2];
  x_Position  = v_Position[2];
  g_clip      = (u_clip == 2) ? v_clip[2] : (u_clip == 1) ? visible( v_Position[2] ) : 0;
  #if WITH_TESSELLATION == 0
  v_Parameter = vec2(0,1);
  #endif

  // i'm not sure if this need to be set after every vertex, or just once before the end of the primitive
  // it at least needs to be set before the last EmitVertex()
  gl_PrimitiveID = gl_PrimitiveIDIn;

  EmitVertex();

  EndPrimitive();

}
