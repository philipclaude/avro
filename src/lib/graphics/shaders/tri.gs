R"(#version 410
layout( triangles ) in;
layout( triangle_strip, max_vertices = 3 ) out;
in vec3 VPosition[];
in vec3 VNormal[];
in vec3 RPosition[];
out vec3 GNormal;
out vec3 GPosition;
out vec3 XYZPosition;
noperspective out vec3 GEdgeDistance;
noperspective out vec3 GTessDistance;
noperspective out vec2 GBarycentric;
uniform vec2  WindowScale;
uniform mat3  NormalMatrix;
uniform mat4  ViewportMatrix;
uniform float shrink;
uniform int   solDeg;
uniform vec4 CutPlane;
uniform vec3  trans;
uniform float scaling;
uniform vec4  scenescaling;
//-- fill in texture array
//out float texture[10];
void main()
{
  // Transform each vertex into viewport space
  //vec2 p0 = vec2(ViewportMatrix * (gl_in[0].gl_Position / gl_in[0].gl_Position.w));
  //vec2 p1 = vec2(ViewportMatrix * (gl_in[1].gl_Position / gl_in[1].gl_Position.w));
  //vec2 p2 = vec2(ViewportMatrix * (gl_in[2].gl_Position / gl_in[2].gl_Position.w));
  //
  //float a = length(p1 - p2);
  //float b = length(p2 - p0);
  //float c = length(p1 - p0);
  //float alpha = acos( (b*b + c*c - a*a) / (2.0*b*c) );
  //float beta = acos( (a*a + c*c - b*b) / (2.0*a*c) );
  //float ha = abs( c * sin( beta ) );
  //float hb = abs( c * sin( alpha ) );
  //float hc = abs( b * sin( alpha ) );
  //float hmax = max(ha,hb);
  //hmax = max(hmax,hc);
  vec4 G  =  0.333333* ( gl_in[0].gl_Position + gl_in[1].gl_Position + gl_in[2].gl_Position );
  vec3 RG =  0.333333* ( RPosition[0] + RPosition[1] + RPosition[2] );
  vec4 vertex_0 = G + shrink*(gl_in[0].gl_Position-G);
  vec4 vertex_1 = G + shrink*(gl_in[1].gl_Position-G);
  vec4 vertex_2 = G + shrink*(gl_in[2].gl_Position-G);
  vec3 A = VPosition[1]-VPosition[0];
  vec3 B = VPosition[2]-VPosition[0];
  vec3 gFacetNormal =  normalize(cross(A, B));
  //vec2 p0 = WindowScale * vertex_0.xy/vertex_0.w;
  //vec2 p1 = WindowScale * vertex_1.xy/vertex_1.w;
  //vec2 p2 = WindowScale * vertex_2.xy/vertex_2.w;
  //
  //vec2 v0 = p2-p1;
  //vec2 v1 = p2-p0;
  //vec2 v2 = p1-p0;
  //float area = abs(v1.x*v2.y - v1.y * v2.x);
  gl_PrimitiveID     = gl_PrimitiveIDIn;
  GPosition          = VPosition[0];
  XYZPosition        = RG + shrink*(RPosition[0] -RG);
  GNormal            = gFacetNormal;
  //GNormal          =   VNormal[0];
  GEdgeDistance      = vec3(1,0,0);
  GTessDistance      = vec3(1,0,0);
  GBarycentric       = vec2(1,0);
  gl_Position        = vertex_0;
	gl_ClipDistance[0] = dot(XYZPosition/scenescaling.w + scenescaling.xyz,CutPlane.xyz) - CutPlane.w;
  EmitVertex();
  gl_PrimitiveID = gl_PrimitiveIDIn;
  GPosition      = VPosition[1];
  XYZPosition    = RG + shrink*(RPosition[1] -RG);
  GNormal        = gFacetNormal;
  //GNormal        =   VNormal[1];
  GEdgeDistance  = vec3(0,1,0);
  GTessDistance  = vec3(0,1,0);
  GBarycentric   = vec2(0,1);
  gl_Position    = vertex_1;
	gl_ClipDistance[0] = dot(XYZPosition/scenescaling.w + scenescaling.xyz,CutPlane.xyz) - CutPlane.w;
  EmitVertex();
  gl_PrimitiveID = gl_PrimitiveIDIn;
  GPosition      = VPosition[2];
  XYZPosition    = RG + shrink*(RPosition[2] -RG);
  GNormal        = gFacetNormal;
  //GNormal        =   VNormal[2];
  GEdgeDistance  = vec3(0,0,1);
  GTessDistance  = vec3(0,0,1);
  GBarycentric   = vec2(0,0);
  gl_Position    = vertex_2;
	gl_ClipDistance[0] = dot(XYZPosition/scenescaling.w + scenescaling.xyz,CutPlane.xyz) - CutPlane.w;
  EmitVertex();
  //int idx    = gl_PrimitiveID*(-(SolDeg - 1)*(SolDeg- 2)*(SolDeg - 3)/6 + 3*SolDeg*(SolDeg- 2)*(SolDeg - 3)/2 - 3*(SolDeg - 1)*SolDeg*(SolDeg - 3) + 5*solDeg*(solDeg -1 )*(solDeg - 2)/3);
  //texture[0] = texelFetch(tex, idx     ).x;
  EndPrimitive();
}
)"
