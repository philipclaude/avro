// the version should be specified by a macro at compile time
// e.g. #version 410 or #version 330

layout (points) in;

uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_NormalMatrix;
uniform mat4 u_ModelViewMatrix;
float u_radius = 0.01; // make uniform later

in vec3 p_Position[];
in float v_Color[];

//out vec3 v_Position;
//out vec3 v_Normal;

out float g_Color;


// total number of points = NU*NV*2*3
#define NU 6
#define NV 5

layout (triangle_strip, max_vertices = 180 ) out;

#define M_PI 3.1415926535897932384626433832795

vec3
spherical( float theta , float phi ) {
  return u_radius * vec3( cos(theta)*sin(phi) , sin(theta)*sin(phi) , cos(phi) );
}

void
main() {

  // index of the particle
  int idx = gl_PrimitiveIDIn;

  int ntheta = NU;
  int nphi   = NV;

  float dtheta = 2.0*M_PI/NU;
  float dphi   = M_PI/NV;

  vec3 normal; // TODO compute sphere normal

  vec3 center = p_Position[0];

  float u0,u1,v0,v1;

  for (int i = 0; i < ntheta; i++) {

    u0 = float(i)*dtheta;
    u1 = float(i+1)*dtheta;

    for (int j = 0; j < nphi; j++) {

      v0 = float(j)*dphi;
      v1 = float(j+1)*dphi;

      // compute coordinates of the vertices
      vec3 p0 = spherical( u0 , v0 );
      vec3 p1 = spherical( u1 , v0 );
      vec3 p2 = spherical( u1 , v1 );
      vec3 p3 = spherical( u0 , v1 );

      gl_Position = u_ModelViewProjectionMatrix * vec4(center + p0 , 1.0 );
      g_Color = v_Color[0];
      EmitVertex();

      gl_Position = u_ModelViewProjectionMatrix * vec4(center + p1 , 1.0 );
      g_Color = v_Color[0];
      EmitVertex();

      gl_Position = u_ModelViewProjectionMatrix * vec4(center + p2 , 1.0 );
      g_Color = v_Color[0];
      gl_PrimitiveID = gl_PrimitiveIDIn;
      EmitVertex();

      EndPrimitive();

      gl_Position = u_ModelViewProjectionMatrix * vec4(center + p0 , 1.0 );
      g_Color = v_Color[0];
      EmitVertex();

      gl_Position = u_ModelViewProjectionMatrix * vec4(center + p2 , 1.0 );
      g_Color = v_Color[0];
      EmitVertex();

      gl_Position = u_ModelViewProjectionMatrix * vec4(center + p3 , 1.0 );
      g_Color = v_Color[0];
      gl_PrimitiveID = gl_PrimitiveIDIn;
      EmitVertex();

      EndPrimitive();

    }
  }

}
