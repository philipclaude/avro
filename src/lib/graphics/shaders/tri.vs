R"(#version 400

layout (location = 0 ) in vec3 VertexPosition;
layout (location = 1 ) in vec3 VertexNormal;

uniform mat3 NormalMatrix;
uniform mat4 ModelViewMatrix;
uniform mat4 MVP;

uniform vec3  trans;
uniform float scaling;

out vec3 VPosition;
out vec3 VNormal;
out vec3 RPosition;

void main()
{
   VNormal         =  NormalMatrix*VertexNormal;
   VPosition       =  vec3(ModelViewMatrix * vec4((VertexPosition),1.0));

   RPosition   = scaling*VertexPosition + trans;
   gl_Position = MVP * vec4(scaling*VertexPosition + trans,1.0);
}
)"
