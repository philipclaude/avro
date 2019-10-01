R"(#version 400

layout( location = 0 ) out vec4 FragColor;

struct LightInfo {
  vec4 Position;  // Light position in eye coords.
  vec3 Intensity; // A,D,S intensity
};

uniform LightInfo Light;

struct MaterialInfo {
  vec3 Ka;            // Ambient reflectivity
  vec3 Kd;            // Diffuse reflectivity
  vec3 Ks;            // Specular reflectivity
  float Shininess;    // Specular shininess factor
};


/* texture access */
uniform samplerBuffer tex;

uniform MaterialInfo Material;

uniform struct LineInfo {
  float Width;
  vec4 Color;
} Line;

in vec3 GPosition;
in vec3 GNormal;
in vec3 XYZPosition;
noperspective in vec3 GEdgeDistance;
noperspective in vec3 GTessDistance;
noperspective in vec2 GBarycentric;

uniform int  NbrHighLight;
uniform int  IdHighlight[10];
uniform bool HighLightAll;
uniform float IsoWidth;

/* display option */
uniform int  BasePicking;
uniform bool DiscardUnPicked;
uniform bool LineOn;
uniform bool TessOn;
uniform bool SolOn;
uniform bool IsoOn;
uniform int  nbrIso;
uniform bool CheckInOut;
uniform bool Picking;
uniform int  UsePhong;
uniform int  TypePrim;
uniform bool FaceOff;
uniform vec4 scenescaling;
uniform float Palette[5];
uniform int   solDeg;
uniform int   Sol;

/*-------------------------------------*/
/*  picking shading: convert id-color  */
/*-------------------------------------*/
vec4 encode_id (int id) {

  int a = id / 16777216;
  int r = (id - a*16777216)/65536;
  int g = (id - a*16777216 - r * 65536) / 256;
  int b = (id - a*16777216 - r * 65536 - g * 256);

  float fr = r;
  float fg = g;
  float fb = b;
  float fa = a;

  vec4 color = vec4( fr / 255.0, fg / 255.0, fb / 255.0, a / 255.0);

  // convert to floats. only divide by 255, because range is 0-255
  return color;
}

/*-------------------------*/
/*  color map function    */
/*-------------------------*/

uniform samplerBuffer colmap;
uniform int colmapsiz;

// --- Return colormap rgb between 0 and 1.
vec3 colormap_rgb(float value, float minVal, float maxVal) {
  vec3 rgb;
  float clamped = clamp(value, minVal, maxVal);
  float step    = (maxVal - minVal) / float(colmapsiz - 1);

  float pos =  fract((clamped - minVal) / step);

  int point1 = int(floor((clamped - minVal) / step));
  int point2 = min( point1 + 1 , colmapsiz - 1);

  rgb.x = mix(texelFetch(colmap, 3*point1).x,     texelFetch(colmap, 3*point2).x,     pos);
  rgb.y = mix(texelFetch(colmap, 3*point1 + 1).x, texelFetch(colmap, 3*point2 + 1).x, pos);
  rgb.z = mix(texelFetch(colmap, 3*point1 + 2).x, texelFetch(colmap, 3*point2 + 2).x, pos);
  return(rgb);
}

/*--------------------------------*/
/*  toon shading (just for fun !) */
/*--------------------------------*/
vec3 LightPosition = vec3(0.25, 0.25, 1);
float Shininess    = 50;

float stepmix(float edge0, float edge1, float E, float x)
{
    float T = clamp(0.5 * (x - edge0 + E) / E, 0.0, 1.0);
    return mix(edge0, edge1, T);
}
vec3 toonModel(vec3 norm)
{
    vec3 N = normalize(norm);
    vec3 L = normalize(LightPosition);
    vec3 Eye = vec3(0, 0, 1);
    vec3 H = normalize(L + Eye);


    float df = max(0.0, dot(N, L));
    float sf = max(0.0, dot(N, H));
    sf = pow(sf, Material.Shininess);

    const float A = 0.1;
    const float B = 0.3;
    const float C = 0.6;
    const float D = 1.0;
    float E = fwidth(df);

    if      (df > A - E && df < A + E) df = stepmix(A, B, E, df);
    else if (df > B - E && df < B + E) df = stepmix(B, C, E, df);
    else if (df > C - E && df < C + E) df = stepmix(C, D, E, df);
    else if (df < A) df = 0.0;
    else if (df < B) df = B;
    else if (df < C) df = C;
    else df = D;

    E = fwidth(sf);
    if (sf > 0.5 - E && sf < 0.5 + E)
    {
        sf = smoothstep(0.5 - E, 0.5 + E, sf);
    }
    else
    {
        sf = step(0.5, sf);
    }

    return Material.Ka + df * Material.Kd + sf * Material.Ks;

}

vec3 toonModelSol(vec3 norm, vec3 color)
{
    vec3 N = normalize(norm);
    vec3 L = normalize(LightPosition);
    vec3 Eye = vec3(0, 0, 1);
    vec3 H = normalize(L + Eye);


    float df = max(0.0, dot(N, L));
    float sf = max(0.0, dot(N, H));
    sf = pow(sf, Material.Shininess);

    const float A = 0.1;
    const float B = 0.3;
    const float C = 0.6;
    const float D = 1.0;
    float E = fwidth(df);

    if      (df > A - E && df < A + E) df = stepmix(A, B, E, df);
    else if (df > B - E && df < B + E) df = stepmix(B, C, E, df);
    else if (df > C - E && df < C + E) df = stepmix(C, D, E, df);
    else if (df < A) df = 0.0;
    else if (df < B) df = B;
    else if (df < C) df = C;
    else df = D;

    E = fwidth(sf);
    if (sf > 0.5 - E && sf < 0.5 + E)
    {
        sf = smoothstep(0.5 - E, 0.5 + E, sf);
    }
    else
    {
        sf = step(0.5, sf);
    }

    return Material.Ka + df * color + sf * Material.Ks;

}

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 rgb2hsv(vec3 c)
{
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

/*--------------------------------*/
/*  diffuse shading (medit-like)  */
/*--------------------------------*/
vec3 diffuseModel( vec3 pos, vec3 norm )
{
    vec3 s    = normalize(vec3(Light.Position) - pos);
    vec3 t    = normalize(norm);
    vec3 spec = Light.Intensity * Material.Kd * max( dot(s,t), 0.0 ) ;

    return  spec;
}

vec3 diffuseModelSol( vec3 pos, vec3 norm , vec3 color)
{
    vec3 s    = normalize(vec3(Light.Position) - pos);
    vec3 t    = normalize(norm);
    vec3 spec = Light.Intensity *color * max( dot(s,t), 0.0 ) ;

    return  spec;
}

/*--------------------------------*/
/*  standard ADS                  */
/*--------------------------------*/
vec3 phongModel( vec3 pos, vec3 norm )
{
    vec3 s = normalize(vec3(Light.Position) - pos);
    vec3 v = normalize(-pos.xyz);
    vec3 r = reflect( -s, norm );
    vec3 ambient = Light.Intensity * Material.Ka;
    float sDotN = max( dot(s,norm), 0.0 );
    vec3 diffuse = Light.Intensity * Material.Kd * sDotN;
    vec3 spec = vec3(0.0);
    if( sDotN > 0.0 )
        spec = Light.Intensity * Material.Ks *
               pow( max( dot(r,v), 0.0 ), Material.Shininess );

    return ambient + diffuse + spec;

}

vec3 phongModelSol( vec3 pos, vec3 norm, vec3 color )
{
    vec3 s = normalize(vec3(Light.Position) - pos);
    vec3 v = normalize(-pos.xyz);
    vec3 r = reflect( -s, norm );
    vec3 ambient = Light.Intensity * Material.Ka;
    float sDotN = max( dot(s,norm), 0.0 );
    vec3 diffuse = Light.Intensity *  color * sDotN;
    vec3 spec = vec3(0.0);
    if( sDotN > 0.0 )
        spec = Light.Intensity * Material.Ks *
               pow( max( dot(r,v), 0.0 ), Material.Shininess );

    return ambient + diffuse + spec;

}


float solp0()
{
 return texelFetch(tex, gl_PrimitiveID).x;
}


float solp1(float u, float v)
{
 int idx    = gl_PrimitiveID*3;
 float p100 = texelFetch(tex, idx     ).x;
 float p010 = texelFetch(tex, idx + 1 ).x;
 float p001 = texelFetch(tex, idx + 2 ).x;
 return  u*p100  + v*p010 + (1-v-u)*p001;
}


float solp2(float u, float v)
{
 int idx    = gl_PrimitiveID*6;
 float p200 = texelFetch(tex, idx     ).x;
 float p020 = texelFetch(tex, idx + 1 ).x;
 float p002 = texelFetch(tex, idx + 2 ).x;
 float p110 = texelFetch(tex, idx + 3 ).x;
 float p011 = texelFetch(tex, idx + 4 ).x;
 float p101 = texelFetch(tex, idx + 5 ).x;

 return  u*(p200*u  + 2.0*p110*v) + p020*v*v +  (1.-u-v)*( 2.0*(p101*u + p011*v) + p002*(1.-u-v) );
}



float solp3(float u, float v)
{

 int idx    = gl_PrimitiveID*10;
 float p300 = texelFetch(tex, idx     ).x;
 float p030 = texelFetch(tex, idx + 1 ).x;
 float p003 = texelFetch(tex, idx + 2 ).x;
 float p210 = texelFetch(tex, idx + 3 ).x;
 float p120 = texelFetch(tex, idx + 4 ).x;
 float p021 = texelFetch(tex, idx + 5 ).x;
 float p012 = texelFetch(tex, idx + 6 ).x;
 float p102 = texelFetch(tex, idx + 7 ).x;
 float p201 = texelFetch(tex, idx + 8 ).x;
 float p111 = texelFetch(tex, idx + 9 ).x;

 return  (u*u*u*p300  + v*v*v*p030 +  (1.-u-v)*(1.-u-v)*(1.-u-v)*p003
        + 3*u*u*v*p210  + 3*u*v*v*p120  + 3*v*v*(1.-u-v)*p021   + 3*v*(1.-u-v)*(1.-u-v)*p012 + 3*u*(1.-u-v)*(1.-u-v)*p102  + 3*u*u*(1.-u-v)*p201
        +  6*u*v*(1.-u-v)*p111);
}

float solp4(float u, float v)
{

 int idx    = gl_PrimitiveID*15;
 float p400 = texelFetch(tex, idx     ).x;
 float p040 = texelFetch(tex, idx + 1 ).x;
 float p004 = texelFetch(tex, idx + 2 ).x;
 float p310 = texelFetch(tex, idx + 3 ).x;
 float p220 = texelFetch(tex, idx + 4 ).x;
 float p130 = texelFetch(tex, idx + 5 ).x;
 float p031 = texelFetch(tex, idx + 6 ).x;
 float p022 = texelFetch(tex, idx + 7 ).x;
 float p013 = texelFetch(tex, idx + 8 ).x;
 float p103 = texelFetch(tex, idx + 9 ).x;
 float p202 = texelFetch(tex, idx + 10 ).x;
 float p301 = texelFetch(tex, idx + 11 ).x;
 float p211 = texelFetch(tex, idx + 12 ).x;
 float p121 = texelFetch(tex, idx + 13 ).x;
 float p112 = texelFetch(tex, idx + 14 ).x;

 return  (u*u*u*u*p400  + v*v*v*v*p040 +  (1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p004
        + 4*u*u*u*v*p310 + 6*u*u*v*v*p220 + 4*u*v*v*v*p130
        + 4*v*v*v*(1.-u-v)*p031  + 6*v*v*(1.-u-v)*(1.-u-v)*p022 + 4*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p013
        + 4*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*p103  + 6*u*u*(1.-u-v)*(1.-u-v)*p202  + 4*u*u*u*(1.-u-v)*p301
        +  12*u*u*v*(1.-u-v)*p211 +   12*u*v*v*(1.-u-v)*p121 +   12*u*v*(1.-u-v)*(1.-u-v)*p112);
}

float solp5(float u, float v)
{

 int idx    = gl_PrimitiveID*21;
 float p500 = texelFetch(tex, idx     ).x;
 float p050 = texelFetch(tex, idx + 1 ).x;
 float p005 = texelFetch(tex, idx + 2 ).x;
 float p410 = texelFetch(tex, idx + 3 ).x;
 float p320 = texelFetch(tex, idx + 4 ).x;
 float p230 = texelFetch(tex, idx + 5 ).x;
 float p140 = texelFetch(tex, idx + 6 ).x;
 float p041 = texelFetch(tex, idx + 7 ).x;
 float p032 = texelFetch(tex, idx + 8 ).x;
 float p023 = texelFetch(tex, idx + 9 ).x;
 float p014 = texelFetch(tex, idx + 10 ).x;
 float p104 = texelFetch(tex, idx + 11 ).x;
 float p203 = texelFetch(tex, idx + 12 ).x;
 float p302 = texelFetch(tex, idx + 13 ).x;
 float p401 = texelFetch(tex, idx + 14 ).x;
 float p311 = texelFetch(tex, idx + 15 ).x;
 float p131 = texelFetch(tex, idx + 16 ).x;
 float p113 = texelFetch(tex, idx + 17 ).x;
 float p221 = texelFetch(tex, idx + 18 ).x;
 float p122 = texelFetch(tex, idx + 19 ).x;
 float p212 = texelFetch(tex, idx + 20 ).x;

 return  (u*u*u*u*u*p500  + v*v*v*v*v*p050 +  (1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p005
        + 5*u*u*u*u*v*p410 + 10*u*u*u*v*v*p320 + 10*u*u*v*v*v*p230  + 5*u*v*v*v*v*p140
        + 5*v*v*v*v*(1.-u-v)*p041  + 10*v*v*v*(1.-u-v)*(1.-u-v)*p032 + 10*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p023 + 5*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p014
        + 5*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p104  + 10*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*p203  + 10*u*u*u*(1.-u-v)*(1.-u-v)*p302  + 5*u*u*u*u*(1.-u-v)*p401
        +  20*u*u*u*v*(1.-u-v)*p311 +   20*u*v*v*v*(1.-u-v)*p131 + 20*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p113
        +  30*u*u*v*v*(1.-u-v)*p221 +   30*u*v*v*(1.-u-v)*(1.-u-v)*p122 + 30*u*u*v*(1.-u-v)*(1.-u-v)*p212);
}

float solp6(float u, float v)
{
 int idx    = gl_PrimitiveID*28;
 float p600 = texelFetch(tex, idx +0 ).x;
 float p060 = texelFetch(tex, idx +1 ).x;
 float p006 = texelFetch(tex, idx +2 ).x;
 float p510 = texelFetch(tex, idx +3 ).x;
 float p420 = texelFetch(tex, idx +4 ).x;
 float p330 = texelFetch(tex, idx +5 ).x;
 float p240 = texelFetch(tex, idx +6 ).x;
 float p150 = texelFetch(tex, idx +7 ).x;
 float p051 = texelFetch(tex, idx +8 ).x;
 float p042 = texelFetch(tex, idx +9 ).x;
 float p033 = texelFetch(tex, idx +10 ).x;
 float p024 = texelFetch(tex, idx +11 ).x;
 float p015 = texelFetch(tex, idx +12 ).x;
 float p105 = texelFetch(tex, idx +13 ).x;
 float p204 = texelFetch(tex, idx +14 ).x;
 float p303 = texelFetch(tex, idx +15 ).x;
 float p402 = texelFetch(tex, idx +16 ).x;
 float p501 = texelFetch(tex, idx +17 ).x;
 float p411 = texelFetch(tex, idx +18 ).x;
 float p141 = texelFetch(tex, idx +19 ).x;
 float p114 = texelFetch(tex, idx +20 ).x;
 float p321 = texelFetch(tex, idx +21 ).x;
 float p231 = texelFetch(tex, idx +22 ).x;
 float p132 = texelFetch(tex, idx +23 ).x;
 float p123 = texelFetch(tex, idx +24 ).x;
 float p213 = texelFetch(tex, idx +25 ).x;
 float p312 = texelFetch(tex, idx +26 ).x;
 float p222 = texelFetch(tex, idx +27 ).x;
 return(u*u*u*u*u*u*p600 + v*v*v*v*v*v*p060 + (1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p006
         + 6*u*u*u*u*u*v*p510 + 15*u*u*u*u*v*v*p420 + 20*u*u*u*v*v*v*p330
         + 15*u*u*v*v*v*v*p240 + 6*u*v*v*v*v*v*p150 + 6*v*v*v*v*v*(1.-u-v)*p051
         + 15*v*v*v*v*(1.-u-v)*(1.-u-v)*p042 + 20*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p033 + 15*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p024
         + 6*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p015 + 6*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p105 + 15*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p204
         + 20*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*p303 + 15*u*u*u*u*(1.-u-v)*(1.-u-v)*p402 + 6*u*u*u*u*u*(1.-u-v)*p501
         + 30*u*u*u*u*v*(1.-u-v)*p411 + 30*u*v*v*v*v*(1.-u-v)*p141 + 30*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p114
         + 60*u*u*u*v*v*(1.-u-v)*p321 + 60*u*u*v*v*v*(1.-u-v)*p231 + 60*u*v*v*v*(1.-u-v)*(1.-u-v)*p132
         + 60*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p123 + 60*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p213 + 60*u*u*u*v*(1.-u-v)*(1.-u-v)*p312
         + 90*u*u*v*v*(1.-u-v)*(1.-u-v)*p222);
}



float solp7(float u, float v)
{
 int idx    = gl_PrimitiveID*36;
 float p700 = texelFetch(tex, idx +0 ).x;
 float p070 = texelFetch(tex, idx +1 ).x;
 float p007 = texelFetch(tex, idx +2 ).x;
 float p610 = texelFetch(tex, idx +3 ).x;
 float p520 = texelFetch(tex, idx +4 ).x;
 float p430 = texelFetch(tex, idx +5 ).x;
 float p340 = texelFetch(tex, idx +6 ).x;
 float p250 = texelFetch(tex, idx +7 ).x;
 float p160 = texelFetch(tex, idx +8 ).x;
 float p061 = texelFetch(tex, idx +9 ).x;
 float p052 = texelFetch(tex, idx +10 ).x;
 float p043 = texelFetch(tex, idx +11 ).x;
 float p034 = texelFetch(tex, idx +12 ).x;
 float p025 = texelFetch(tex, idx +13 ).x;
 float p016 = texelFetch(tex, idx +14 ).x;
 float p106 = texelFetch(tex, idx +15 ).x;
 float p205 = texelFetch(tex, idx +16 ).x;
 float p304 = texelFetch(tex, idx +17 ).x;
 float p403 = texelFetch(tex, idx +18 ).x;
 float p502 = texelFetch(tex, idx +19 ).x;
 float p601 = texelFetch(tex, idx +20 ).x;
 float p511 = texelFetch(tex, idx +21 ).x;
 float p151 = texelFetch(tex, idx +22 ).x;
 float p115 = texelFetch(tex, idx +23 ).x;
 float p421 = texelFetch(tex, idx +24 ).x;
 float p331 = texelFetch(tex, idx +25 ).x;
 float p241 = texelFetch(tex, idx +26 ).x;
 float p142 = texelFetch(tex, idx +27 ).x;
 float p133 = texelFetch(tex, idx +28 ).x;
 float p124 = texelFetch(tex, idx +29 ).x;
 float p214 = texelFetch(tex, idx +30 ).x;
 float p313 = texelFetch(tex, idx +31 ).x;
 float p412 = texelFetch(tex, idx +32 ).x;
 float p322 = texelFetch(tex, idx +33 ).x;
 float p232 = texelFetch(tex, idx +34 ).x;
 float p223 = texelFetch(tex, idx +35 ).x;
 return(u*u*u*u*u*u*u*p700 + v*v*v*v*v*v*v*p070 + (1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p007
         + 7*u*u*u*u*u*u*v*p610 + 21*u*u*u*u*u*v*v*p520 + 35*u*u*u*u*v*v*v*p430
         + 35*u*u*u*v*v*v*v*p340 + 21*u*u*v*v*v*v*v*p250 + 7*u*v*v*v*v*v*v*p160
         + 7*v*v*v*v*v*v*(1.-u-v)*p061 + 21*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p052 + 35*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p043
         + 35*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p034 + 21*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p025 + 7*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p016
         + 7*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p106 + 21*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p205 + 35*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p304
         + 35*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*p403 + 21*u*u*u*u*u*(1.-u-v)*(1.-u-v)*p502 + 7*u*u*u*u*u*u*(1.-u-v)*p601
         + 42*u*u*u*u*u*v*(1.-u-v)*p511 + 42*u*v*v*v*v*v*(1.-u-v)*p151 + 42*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p115
         + 105*u*u*u*u*v*v*(1.-u-v)*p421 + 140*u*u*u*v*v*v*(1.-u-v)*p331 + 105*u*u*v*v*v*v*(1.-u-v)*p241
         + 105*u*v*v*v*v*(1.-u-v)*(1.-u-v)*p142 + 140*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p133 + 105*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p124
         + 105*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p214 + 140*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p313 + 105*u*u*u*u*v*(1.-u-v)*(1.-u-v)*p412
         + 210*u*u*u*v*v*(1.-u-v)*(1.-u-v)*p322 + 210*u*u*v*v*v*(1.-u-v)*(1.-u-v)*p232 + 210*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p223);
}

float solp8(float u, float v)
{
 int idx    = gl_PrimitiveID*45;
 float p800 = texelFetch(tex, idx +0 ).x;
 float p080 = texelFetch(tex, idx +1 ).x;
 float p008 = texelFetch(tex, idx +2 ).x;
 float p710 = texelFetch(tex, idx +3 ).x;
 float p620 = texelFetch(tex, idx +4 ).x;
 float p530 = texelFetch(tex, idx +5 ).x;
 float p440 = texelFetch(tex, idx +6 ).x;
 float p350 = texelFetch(tex, idx +7 ).x;
 float p260 = texelFetch(tex, idx +8 ).x;
 float p170 = texelFetch(tex, idx +9 ).x;
 float p071 = texelFetch(tex, idx +10 ).x;
 float p062 = texelFetch(tex, idx +11 ).x;
 float p053 = texelFetch(tex, idx +12 ).x;
 float p044 = texelFetch(tex, idx +13 ).x;
 float p035 = texelFetch(tex, idx +14 ).x;
 float p026 = texelFetch(tex, idx +15 ).x;
 float p017 = texelFetch(tex, idx +16 ).x;
 float p107 = texelFetch(tex, idx +17 ).x;
 float p206 = texelFetch(tex, idx +18 ).x;
 float p305 = texelFetch(tex, idx +19 ).x;
 float p404 = texelFetch(tex, idx +20 ).x;
 float p503 = texelFetch(tex, idx +21 ).x;
 float p602 = texelFetch(tex, idx +22 ).x;
 float p701 = texelFetch(tex, idx +23 ).x;
 float p611 = texelFetch(tex, idx +24 ).x;
 float p161 = texelFetch(tex, idx +25 ).x;
 float p116 = texelFetch(tex, idx +26 ).x;
 float p521 = texelFetch(tex, idx +27 ).x;
 float p431 = texelFetch(tex, idx +28 ).x;
 float p341 = texelFetch(tex, idx +29 ).x;
 float p251 = texelFetch(tex, idx +30 ).x;
 float p152 = texelFetch(tex, idx +31 ).x;
 float p143 = texelFetch(tex, idx +32 ).x;
 float p134 = texelFetch(tex, idx +33 ).x;
 float p125 = texelFetch(tex, idx +34 ).x;
 float p215 = texelFetch(tex, idx +35 ).x;
 float p314 = texelFetch(tex, idx +36 ).x;
 float p413 = texelFetch(tex, idx +37 ).x;
 float p512 = texelFetch(tex, idx +38 ).x;
 float p422 = texelFetch(tex, idx +39 ).x;
 float p242 = texelFetch(tex, idx +40 ).x;
 float p224 = texelFetch(tex, idx +41 ).x;
 float p332 = texelFetch(tex, idx +42 ).x;
 float p233 = texelFetch(tex, idx +43 ).x;
 float p323 = texelFetch(tex, idx +44 ).x;
 return(u*u*u*u*u*u*u*u*p800 + v*v*v*v*v*v*v*v*p080 + (1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p008
         + 8*u*u*u*u*u*u*u*v*p710 + 28*u*u*u*u*u*u*v*v*p620 + 56*u*u*u*u*u*v*v*v*p530
         + 70*u*u*u*u*v*v*v*v*p440 + 56*u*u*u*v*v*v*v*v*p350 + 28*u*u*v*v*v*v*v*v*p260
         + 8*u*v*v*v*v*v*v*v*p170 + 8*v*v*v*v*v*v*v*(1.-u-v)*p071 + 28*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p062
         + 56*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p053 + 70*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p044 + 56*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p035
         + 28*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p026 + 8*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p017 + 8*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p107
         + 28*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p206 + 56*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p305 + 70*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p404
         + 56*u*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*p503 + 28*u*u*u*u*u*u*(1.-u-v)*(1.-u-v)*p602 + 8*u*u*u*u*u*u*u*(1.-u-v)*p701
         + 56*u*u*u*u*u*u*v*(1.-u-v)*p611 + 56*u*v*v*v*v*v*v*(1.-u-v)*p161 + 56*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p116
         + 168*u*u*u*u*u*v*v*(1.-u-v)*p521 + 280*u*u*u*u*v*v*v*(1.-u-v)*p431 + 280*u*u*u*v*v*v*v*(1.-u-v)*p341
         + 168*u*u*v*v*v*v*v*(1.-u-v)*p251 + 168*u*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p152 + 280*u*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p143
         + 280*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p134 + 168*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p125 + 168*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p215
         + 280*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p314 + 280*u*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p413 + 168*u*u*u*u*u*v*(1.-u-v)*(1.-u-v)*p512
         + 420*u*u*u*u*v*v*(1.-u-v)*(1.-u-v)*p422 + 420*u*u*v*v*v*v*(1.-u-v)*(1.-u-v)*p242 + 420*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p224
         + 560*u*u*u*v*v*v*(1.-u-v)*(1.-u-v)*p332 + 560*u*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p233 + 560*u*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p323);
}

float solp9(float u, float v)
{
 int idx    = gl_PrimitiveID*55;
 float p900 = texelFetch(tex, idx +0 ).x;
 float p090 = texelFetch(tex, idx +1 ).x;
 float p009 = texelFetch(tex, idx +2 ).x;
 float p810 = texelFetch(tex, idx +3 ).x;
 float p720 = texelFetch(tex, idx +4 ).x;
 float p630 = texelFetch(tex, idx +5 ).x;
 float p540 = texelFetch(tex, idx +6 ).x;
 float p450 = texelFetch(tex, idx +7 ).x;
 float p360 = texelFetch(tex, idx +8 ).x;
 float p270 = texelFetch(tex, idx +9 ).x;
 float p180 = texelFetch(tex, idx +10 ).x;
 float p081 = texelFetch(tex, idx +11 ).x;
 float p072 = texelFetch(tex, idx +12 ).x;
 float p063 = texelFetch(tex, idx +13 ).x;
 float p054 = texelFetch(tex, idx +14 ).x;
 float p045 = texelFetch(tex, idx +15 ).x;
 float p036 = texelFetch(tex, idx +16 ).x;
 float p027 = texelFetch(tex, idx +17 ).x;
 float p018 = texelFetch(tex, idx +18 ).x;
 float p108 = texelFetch(tex, idx +19 ).x;
 float p207 = texelFetch(tex, idx +20 ).x;
 float p306 = texelFetch(tex, idx +21 ).x;
 float p405 = texelFetch(tex, idx +22 ).x;
 float p504 = texelFetch(tex, idx +23 ).x;
 float p603 = texelFetch(tex, idx +24 ).x;
 float p702 = texelFetch(tex, idx +25 ).x;
 float p801 = texelFetch(tex, idx +26 ).x;
 float p711 = texelFetch(tex, idx +27 ).x;
 float p171 = texelFetch(tex, idx +28 ).x;
 float p117 = texelFetch(tex, idx +29 ).x;
 float p621 = texelFetch(tex, idx +30 ).x;
 float p531 = texelFetch(tex, idx +31 ).x;
 float p441 = texelFetch(tex, idx +32 ).x;
 float p351 = texelFetch(tex, idx +33 ).x;
 float p261 = texelFetch(tex, idx +34 ).x;
 float p162 = texelFetch(tex, idx +35 ).x;
 float p153 = texelFetch(tex, idx +36 ).x;
 float p144 = texelFetch(tex, idx +37 ).x;
 float p135 = texelFetch(tex, idx +38 ).x;
 float p126 = texelFetch(tex, idx +39 ).x;
 float p216 = texelFetch(tex, idx +40 ).x;
 float p315 = texelFetch(tex, idx +41 ).x;
 float p414 = texelFetch(tex, idx +42 ).x;
 float p513 = texelFetch(tex, idx +43 ).x;
 float p612 = texelFetch(tex, idx +44 ).x;
 float p522 = texelFetch(tex, idx +45 ).x;
 float p252 = texelFetch(tex, idx +46 ).x;
 float p225 = texelFetch(tex, idx +47 ).x;
 float p432 = texelFetch(tex, idx +48 ).x;
 float p342 = texelFetch(tex, idx +49 ).x;
 float p243 = texelFetch(tex, idx +50 ).x;
 float p234 = texelFetch(tex, idx +51 ).x;
 float p324 = texelFetch(tex, idx +52 ).x;
 float p423 = texelFetch(tex, idx +53 ).x;
 float p333 = texelFetch(tex, idx +54 ).x;
 return(u*u*u*u*u*u*u*u*u*p900 + v*v*v*v*v*v*v*v*v*p090 + (1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p009
         + 9*u*u*u*u*u*u*u*u*v*p810 + 36*u*u*u*u*u*u*u*v*v*p720 + 84*u*u*u*u*u*u*v*v*v*p630
         + 126*u*u*u*u*u*v*v*v*v*p540 + 126*u*u*u*u*v*v*v*v*v*p450 + 84*u*u*u*v*v*v*v*v*v*p360
         + 36*u*u*v*v*v*v*v*v*v*p270 + 9*u*v*v*v*v*v*v*v*v*p180 + 9*v*v*v*v*v*v*v*v*(1.-u-v)*p081
         + 36*v*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p072 + 84*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p063 + 126*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p054
         + 126*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p045 + 84*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p036 + 36*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p027
         + 9*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p018 + 9*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p108 + 36*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p207
         + 84*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p306 + 126*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p405 + 126*u*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p504
         + 84*u*u*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*p603 + 36*u*u*u*u*u*u*u*(1.-u-v)*(1.-u-v)*p702 + 9*u*u*u*u*u*u*u*u*(1.-u-v)*p801
         + 72*u*u*u*u*u*u*u*v*(1.-u-v)*p711 + 72*u*v*v*v*v*v*v*v*(1.-u-v)*p171 + 72*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p117
         + 252*u*u*u*u*u*u*v*v*(1.-u-v)*p621 + 504*u*u*u*u*u*v*v*v*(1.-u-v)*p531 + 630*u*u*u*u*v*v*v*v*(1.-u-v)*p441
         + 504*u*u*u*v*v*v*v*v*(1.-u-v)*p351 + 252*u*u*v*v*v*v*v*v*(1.-u-v)*p261 + 252*u*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p162
         + 504*u*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p153 + 630*u*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p144 + 504*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p135
         + 252*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p126 + 252*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p216 + 504*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p315
         + 630*u*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p414 + 504*u*u*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p513 + 252*u*u*u*u*u*u*v*(1.-u-v)*(1.-u-v)*p612
         + 756*u*u*u*u*u*v*v*(1.-u-v)*(1.-u-v)*p522 + 756*u*u*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p252 + 756*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p225
         + 1260*u*u*u*u*v*v*v*(1.-u-v)*(1.-u-v)*p432 + 1260*u*u*u*v*v*v*v*(1.-u-v)*(1.-u-v)*p342 + 1260*u*u*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p243
         + 1260*u*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p234 + 1260*u*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p324 + 1260*u*u*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p423
         + 1680*u*u*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p333);
}



float solp10(float u, float v)
{
 int idx    = gl_PrimitiveID*66;
 float p1000 = texelFetch(tex, idx +0 ).x;
 float p0100 = texelFetch(tex, idx +1 ).x;
 float p0010 = texelFetch(tex, idx +2 ).x;
 float p910 = texelFetch(tex, idx +3 ).x;
 float p820 = texelFetch(tex, idx +4 ).x;
 float p730 = texelFetch(tex, idx +5 ).x;
 float p640 = texelFetch(tex, idx +6 ).x;
 float p550 = texelFetch(tex, idx +7 ).x;
 float p460 = texelFetch(tex, idx +8 ).x;
 float p370 = texelFetch(tex, idx +9 ).x;
 float p280 = texelFetch(tex, idx +10 ).x;
 float p190 = texelFetch(tex, idx +11 ).x;
 float p091 = texelFetch(tex, idx +12 ).x;
 float p082 = texelFetch(tex, idx +13 ).x;
 float p073 = texelFetch(tex, idx +14 ).x;
 float p064 = texelFetch(tex, idx +15 ).x;
 float p055 = texelFetch(tex, idx +16 ).x;
 float p046 = texelFetch(tex, idx +17 ).x;
 float p037 = texelFetch(tex, idx +18 ).x;
 float p028 = texelFetch(tex, idx +19 ).x;
 float p019 = texelFetch(tex, idx +20 ).x;
 float p109 = texelFetch(tex, idx +21 ).x;
 float p208 = texelFetch(tex, idx +22 ).x;
 float p307 = texelFetch(tex, idx +23 ).x;
 float p406 = texelFetch(tex, idx +24 ).x;
 float p505 = texelFetch(tex, idx +25 ).x;
 float p604 = texelFetch(tex, idx +26 ).x;
 float p703 = texelFetch(tex, idx +27 ).x;
 float p802 = texelFetch(tex, idx +28 ).x;
 float p901 = texelFetch(tex, idx +29 ).x;
 float p811 = texelFetch(tex, idx +30 ).x;
 float p181 = texelFetch(tex, idx +31 ).x;
 float p118 = texelFetch(tex, idx +32 ).x;
 float p721 = texelFetch(tex, idx +33 ).x;
 float p631 = texelFetch(tex, idx +34 ).x;
 float p541 = texelFetch(tex, idx +35 ).x;
 float p451 = texelFetch(tex, idx +36 ).x;
 float p361 = texelFetch(tex, idx +37 ).x;
 float p271 = texelFetch(tex, idx +38 ).x;
 float p172 = texelFetch(tex, idx +39 ).x;
 float p163 = texelFetch(tex, idx +40 ).x;
 float p154 = texelFetch(tex, idx +41 ).x;
 float p145 = texelFetch(tex, idx +42 ).x;
 float p136 = texelFetch(tex, idx +43 ).x;
 float p127 = texelFetch(tex, idx +44 ).x;
 float p217 = texelFetch(tex, idx +45 ).x;
 float p316 = texelFetch(tex, idx +46 ).x;
 float p415 = texelFetch(tex, idx +47 ).x;
 float p514 = texelFetch(tex, idx +48 ).x;
 float p613 = texelFetch(tex, idx +49 ).x;
 float p712 = texelFetch(tex, idx +50 ).x;
 float p622 = texelFetch(tex, idx +51 ).x;
 float p262 = texelFetch(tex, idx +52 ).x;
 float p226 = texelFetch(tex, idx +53 ).x;
 float p532 = texelFetch(tex, idx +54 ).x;
 float p442 = texelFetch(tex, idx +55 ).x;
 float p352 = texelFetch(tex, idx +56 ).x;
 float p253 = texelFetch(tex, idx +57 ).x;
 float p244 = texelFetch(tex, idx +58 ).x;
 float p235 = texelFetch(tex, idx +59 ).x;
 float p325 = texelFetch(tex, idx +60 ).x;
 float p424 = texelFetch(tex, idx +61 ).x;
 float p523 = texelFetch(tex, idx +62 ).x;
 float p433 = texelFetch(tex, idx +63 ).x;
 float p343 = texelFetch(tex, idx +64 ).x;
 float p334 = texelFetch(tex, idx +65 ).x;
 return(u*u*u*u*u*u*u*u*u*u*p1000 + v*v*v*v*v*v*v*v*v*v*p0100 + (1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p0010
         + 10*u*u*u*u*u*u*u*u*u*v*p910 + 45*u*u*u*u*u*u*u*u*v*v*p820 + 120*u*u*u*u*u*u*u*v*v*v*p730
         + 210*u*u*u*u*u*u*v*v*v*v*p640 + 252*u*u*u*u*u*v*v*v*v*v*p550 + 210*u*u*u*u*v*v*v*v*v*v*p460
         + 120*u*u*u*v*v*v*v*v*v*v*p370 + 45*u*u*v*v*v*v*v*v*v*v*p280 + 10*u*v*v*v*v*v*v*v*v*v*p190
         + 10*v*v*v*v*v*v*v*v*v*(1.-u-v)*p091 + 45*v*v*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p082 + 120*v*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p073
         + 210*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p064 + 252*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p055 + 210*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p046
         + 120*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p037 + 45*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p028 + 10*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p019
         + 10*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p109 + 45*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p208 + 120*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p307
         + 210*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p406 + 252*u*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p505 + 210*u*u*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p604
         + 120*u*u*u*u*u*u*u*(1.-u-v)*(1.-u-v)*(1.-u-v)*p703 + 45*u*u*u*u*u*u*u*u*(1.-u-v)*(1.-u-v)*p802 + 10*u*u*u*u*u*u*u*u*u*(1.-u-v)*p901
         + 90*u*u*u*u*u*u*u*u*v*(1.-u-v)*p811 + 90*u*v*v*v*v*v*v*v*v*(1.-u-v)*p181 + 90*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p118
         + 360*u*u*u*u*u*u*u*v*v*(1.-u-v)*p721 + 840*u*u*u*u*u*u*v*v*v*(1.-u-v)*p631 + 1260*u*u*u*u*u*v*v*v*v*(1.-u-v)*p541
         + 1260*u*u*u*u*v*v*v*v*v*(1.-u-v)*p451 + 840*u*u*u*v*v*v*v*v*v*(1.-u-v)*p361 + 360*u*u*v*v*v*v*v*v*v*(1.-u-v)*p271
         + 360*u*v*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p172 + 840*u*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p163 + 1260*u*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p154
         + 1260*u*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p145 + 840*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p136 + 360*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p127
         + 360*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p217 + 840*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p316 + 1260*u*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p415
         + 1260*u*u*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p514 + 840*u*u*u*u*u*u*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p613 + 360*u*u*u*u*u*u*u*v*(1.-u-v)*(1.-u-v)*p712
         + 1260*u*u*u*u*u*u*v*v*(1.-u-v)*(1.-u-v)*p622 + 1260*u*u*v*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p262 + 1260*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p226
         + 2520*u*u*u*u*u*v*v*v*(1.-u-v)*(1.-u-v)*p532 + 3150*u*u*u*u*v*v*v*v*(1.-u-v)*(1.-u-v)*p442 + 2520*u*u*u*v*v*v*v*v*(1.-u-v)*(1.-u-v)*p352
         + 2520*u*u*v*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p253 + 3150*u*u*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p244 + 2520*u*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p235
         + 2520*u*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p325 + 3150*u*u*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p424 + 2520*u*u*u*u*u*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p523
         + 4200*u*u*u*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p433 + 4200*u*u*u*v*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*p343 + 4200*u*u*u*v*v*v*(1.-u-v)*(1.-u-v)*(1.-u-v)*(1.-u-v)*p334);

  //-- Horner form
  // return(1*p0010+(-10*p0010+10*p019+(45*p0010-90*p019+45*p028+(-120*p0010+360*p019-360*p028+120*p037+(210*p0010-840*p019+1260*p028+210*p046-840*p037+
  //      (1260*p019+2520*p037-2520*p028-252*p0010+252*p055-1260*p046+(-1260*p055-4200*p037+3150*p028+3150*p046+210*p0010-1260*p019+210*p064+(4200*p037+2520*p055+
  //      120*p073-2520*p028+840*p019-840*p064-120*p0010-4200*p046+(45*p082+1260*p028-2520*p037+45*p0010-2520*p055-360*p019+3150*p046-360*p073+1260*p064+
  //      (1260*p055+10*p091-360*p028-10*p0010-90*p082+360*p073+840*p037-1260*p046+90*p019-840*p064+(45*p082+p0100+210*p046+1*p0010-10*p091+45*p028-120*p037
  //      -10*p019-120*p073-252*p055+210*p064)*v)*v)*v)*v)*v)*v)*v)*v)*v)*v+(-10*p0010+10*p109+(-90*p109+90*p118-90*p019+90*p0010+(-360*p028+360*p127-360*p0010+
  //      360*p109-720*p118+720*p019+(840*p136+840*p0010-840*p109+2520*p118-2520*p019+2520*p028-840*p037-2520*p127+(1260*p109-1260*p046-7560*p028+5040*p037+
  //      1260*p145+7560*p127-5040*p136-1260*p0010+5040*p019-5040*p118+(12600*p028-6300*p019+6300*p046-12600*p037-6300*p145-12600*p127+6300*p118+1260*p0010
  //      -1260*p055+1260*p154+12600*p136-1260*p109+(-12600*p028+12600*p127+840*p163-840*p0010+5040*p055-840*p064-5040*p154-16800*p136+5040*p019+12600*p145+
  //      840*p109+16800*p037-5040*p118-12600*p046+(-2520*p163+2520*p064+360*p172-360*p073-12600*p145-7560*p127-7560*p055-360*p109+7560*p028+360*p0010+
  //      12600*p136+12600*p046-2520*p019+7560*p154-12600*p037+2520*p118+(6300*p145+720*p019-2520*p028+2520*p163-720*p172-2520*p064-5040*p154-5040*p136+
  //      90*p109-90*p0010+720*p073+90*p181-6300*p046-720*p118-90*p082+2520*p127+5040*p037+5040*p055+(-840*p163+1260*p154+840*p136+1260*p046-1260*p145+
  //      10*p0010-360*p073-840*p037-90*p019+360*p028+840*p064-90*p181+90*p118+10*p190+90*p082-360*p127-10*p091-10*p109-1260*p055+360*p172)*v)*v)*v)*v)*v)*v)*v)*v)*v
  //      +(-90*p109+45*p208+45*p0010+(720*p109-720*p118+360*p019-360*p208+360*p217-360*p0010+(5040*p118+1260*p0010+1260*p208+1260*p226+1260*p028-2520*p109
  //      -2520*p127-2520*p217-2520*p019+(-7560*p226-15120*p118-2520*p0010+5040*p109-5040*p136+15120*p127+2520*p037+7560*p019+2520*p235-2520*p208+7560*p217
  //      -7560*p028+(3150*p0010+18900*p226-12600*p235+3150*p244+25200*p136+25200*p118-37800*p127-12600*p217+18900*p028+3150*p046+3150*p208-12600*p019-6300*p109
  //      -6300*p145-12600*p037+(50400*p127-2520*p0010-12600*p046-25200*p226-12600*p244+2520*p253+2520*p055+25200*p235+25200*p145-5040*p154-2520*p208-25200*p118
  //      -50400*p136-25200*p028+12600*p019+12600*p217+5040*p109+25200*p037+(-2520*p163+1260*p0010+18900*p226-25200*p235+18900*p244-37800*p127+1260*p064+15120*p118
  //      -37800*p145+1260*p208-7560*p217+50400*p136+18900*p046-7560*p253-7560*p019+1260*p262-2520*p109+18900*p028+15120*p154-7560*p055-25200*p037+(-12600*p244
  //      +7560*p253+2520*p217-360*p208+720*p109+2520*p019-7560*p226-2520*p262-7560*p028+12600*p037-12600*p046+7560*p055-2520*p064+360*p073+12600*p235-360*p0010
  //      -5040*p118+5040*p163-720*p172+360*p271+(3150*p244-2520*p253-360*p217+45*p208-90*p109-360*p019+1260*p226+1260*p262+1260*p028-2520*p037+3150*p046
  //      -2520*p055+1260*p064-360*p073+45*p082-2520*p235+45*p280+45*p0010+720*p118-90*p181-2520*p163+720*p172-360*p271+5040*p154-6300*p145+5040*p136-2520*p127)*v
  //      -15120*p154+25200*p145-25200*p136+15120*p127)*v)*v)*v)*v)*v)*v)*v+(120*p307-360*p208-120*p0010+360*p109+(-840*p307-2520*p217+840*p0010+840*p316+2520*p208
  //      -840*p019-2520*p109+2520*p118+(2520*p307-5040*p316+5040*p019-2520*p0010+7560*p127+7560*p109-15120*p118-2520*p028+15120*p217-7560*p226+2520*p325-7560*p208
  //      +(-12600*p235+4200*p334-12600*p109-4200*p037+12600*p028+12600*p316-4200*p307-12600*p019-12600*p325+4200*p0010+12600*p136+37800*p226-37800*p127-37800*p217
  //      +37800*p118+12600*p208+(-75600*p226+12600*p109+50400*p235+16800*p037+12600*p145+4200*p343+16800*p019-50400*p118+4200*p307-12600*p208-25200*p028-50400*p136
  //      -16800*p316+75600*p127-16800*p334-4200*p0010-12600*p244-4200*p046+50400*p217+25200*p325+(2520*p352+37800*p244-7560*p253+12600*p316-37800*p217-2520*p307
  //      +7560*p208-7560*p109-12600*p019+75600*p226+25200*p028-25200*p037+12600*p046-2520*p055-75600*p235+2520*p0010+37800*p118-12600*p343-25200*p325+(-5040*p352
  //      -37800*p244+15120*p253-5040*p316+15120*p217+840*p307-2520*p208+2520*p109+5040*p019-37800*p226-2520*p262-12600*p028+16800*p037-12600*p046+5040*p055
  //      -840*p064+50400*p235-840*p0010-15120*p118+2520*p163+840*p361+12600*p343+12600*p325+(2520*p352+12600*p244-7560*p253+840*p316-2520*p217-120*p307+360*p208
  //      -360*p109-840*p019+7560*p226+2520*p262+2520*p028-4200*p037+4200*p046-2520*p055+840*p064-120*p073-12600*p235+120*p370+120*p0010+2520*p118-2520*p163
  //      +360*p172-360*p271-840*p361-4200*p343-2520*p325+4200*p334+7560*p154-12600*p145+12600*p136-7560*p127)*v-16800*p334-15120*p154+37800*p145-50400*p136
  //      +37800*p127)*v+25200*p334+7560*p154-37800*p145+75600*p136-75600*p127)*v)*v)*v)*v)*v+(-840*p109+210*p0010-840*p307+1260*p208+210*p406+(-7560*p208-5040*p118
  //      +5040*p109+1260*p019-5040*p316-1260*p0010+5040*p307+7560*p217+1260*p415-1260*p406+(-12600*p307+18900*p208+3150*p424+3150*p0010+25200*p316-12600*p325
  //      -12600*p109+18900*p226+25200*p118+3150*p028-6300*p415-37800*p217-12600*p127+3150*p406-6300*p019+(16800*p307-12600*p424+16800*p109+25200*p235+50400*p325
  //      +4200*p037+12600*p415+4200*p433-12600*p028+50400*p127-50400*p316-16800*p334-50400*p118+75600*p217-4200*p406+12600*p019-25200*p208-75600*p226-4200*p0010
  //      -16800*p136+(3150*p442+18900*p244-12600*p415+50400*p316-75600*p217-12600*p433-12600*p307+18900*p208-12600*p109-12600*p019+18900*p028-12600*p037+3150*p046
  //      -75600*p235+3150*p0010+50400*p118+18900*p424-12600*p343-75600*p325+3150*p406+(-5040*p352-6300*p442-37800*p244+6300*p415+7560*p253-25200*p316+37800*p217
  //      +12600*p433+5040*p307-7560*p208+5040*p109+6300*p019-75600*p226-12600*p028+12600*p037-6300*p046+1260*p055+75600*p235-1260*p0010-25200*p118+1260*p451
  //      -12600*p424+25200*p343+50400*p325-1260*p406+(5040*p352+3150*p442+18900*p244-1260*p415-7560*p253+5040*p316-7560*p217-4200*p433-840*p307+1260*p208-840*p109
  //      -1260*p019+18900*p226+1260*p262+3150*p028-4200*p037+3150*p046-1260*p055+210*p064-25200*p235+210*p460+210*p0010+5040*p118-840*p163-840*p361-1260*p451
  //      +3150*p424-12600*p343-12600*p325+210*p406+16800*p334+5040*p154-12600*p145+16800*p136-12600*p127)*v-50400*p334-5040*p154+25200*p145-50400*p136+50400*p127)*v
  //      +50400*p334-12600*p145+50400*p136-75600*p127+113400*p226)*v)*v)*v)*v+(-2520*p208+2520*p307+1260*p109+252*p505-1260*p406-252*p0010+(12600*p208+1260*p514
  //      -6300*p415+6300*p406+12600*p316-12600*p217-6300*p109+1260*p0010-12600*p307+6300*p118-1260*p505-1260*p019+(25200*p415+2520*p505-5040*p514+25200*p307+12600*p109
  //      -25200*p208-50400*p316+2520*p523+50400*p217-2520*p0010-25200*p118-25200*p226-2520*p028+25200*p325-12600*p406+5040*p019-12600*p424+12600*p127+(7560*p514
  //      -37800*p415+75600*p316-75600*p217-12600*p433-25200*p307+25200*p208-12600*p109-7560*p019+2520*p532+75600*p226+7560*p028-2520*p037-25200*p235+2520*p0010
  //      -7560*p523+37800*p118+37800*p424-75600*p325-2520*p505+12600*p406+(-6300*p442-12600*p244-5040*p514+25200*p415-50400*p316+50400*p217+25200*p433+12600*p307
  //      -12600*p208+6300*p109+5040*p019-5040*p532-75600*p226-7560*p028+5040*p037-1260*p046+50400*p235-1260*p0010+7560*p523-25200*p118+1260*p541-37800*p424
  //      +12600*p343+75600*p325+1260*p505-6300*p406+(2520*p352+6300*p442+12600*p244+1260*p514-6300*p415-2520*p253+12600*p316-12600*p217-12600*p433-2520*p307
  //      +2520*p208-1260*p109-1260*p019+2520*p532+25200*p226+2520*p028-2520*p037+1260*p046-252*p055-25200*p235+252*p550+252*p0010-2520*p523+6300*p118-1260*p451
  //      -1260*p541+12600*p424-12600*p343-25200*p325-252*p505+1260*p406+25200*p334+1260*p154-6300*p145+12600*p136-12600*p127)*v-50400*p334+6300*p145-25200*p136
  //      +37800*p127)*v+25200*p334+12600*p136-37800*p127)*v)*v)*v+(-1260*p505-4200*p307-1260*p109+3150*p406+210*p604+210*p0010+3150*p208+(5040*p505-12600*p406
  //      +840*p019-16800*p316+12600*p415-12600*p208+12600*p217+16800*p307-840*p0010+5040*p109-5040*p118-5040*p514+840*p613-840*p604+(18900*p208-37800*p415-25200*p307
  //      +50400*p316+18900*p424-37800*p217+1260*p604+1260*p0010-7560*p505-7560*p127+18900*p406+18900*p226-2520*p613-7560*p109-7560*p523-25200*p325-2520*p019+15120*p118
  //      +15120*p514+1260*p028+1260*p622+(-15120*p514+37800*p415-50400*p316+37800*p217+12600*p433+16800*p307-12600*p208+5040*p109+2520*p019-5040*p532-37800*p226
  //      -2520*p622-2520*p028+840*p037+12600*p235-840*p0010+15120*p523+840*p631-15120*p118-37800*p424+2520*p613+50400*p325-840*p604+5040*p505-12600*p406
  //      +(3150*p442+3150*p244+5040*p514-12600*p415+16800*p316-12600*p217-12600*p433-4200*p307+3150*p208-1260*p109-840*p019+5040*p532+18900*p226+1260*p622
  //      +1260*p028-840*p037+210*p046-12600*p235+210*p640+210*p0010-7560*p523-840*p631+5040*p118-1260*p541+18900*p424-4200*p343-840*p613-25200*p325+210*p604
  //      -1260*p505+3150*p406+16800*p334-1260*p145+5040*p136-7560*p127)*v-16800*p334-5040*p136+15120*p127)*v)*v)*v+(4200*p307-4200*p406-120*p0010+840*p109-840*p604
  //      +120*p703-2520*p208+2520*p505+(360*p712+360*p0010-2520*p613-12600*p415+12600*p406-360*p703+12600*p316+2520*p118+2520*p604-360*p019+7560*p208-7560*p505
  //      -12600*p307-2520*p109-7560*p217+7560*p514+(-15120*p514+25200*p415-25200*p316+15120*p217+12600*p307-7560*p208+2520*p109+720*p019-7560*p226-2520*p622-720*p712
  //      -360*p028-360*p0010+7560*p523+360*p721-5040*p118-12600*p424+5040*p613+360*p703+12600*p325-2520*p604+7560*p505-12600*p406+(7560*p514-12600*p415+12600*p316
  //      -7560*p217-4200*p433-4200*p307+2520*p208-840*p109-360*p019+2520*p532+7560*p226+2520*p622+360*p712+360*p028-120*p037-2520*p235+120*p730+120*p0010-7560*p523
  //      -840*p631-360*p721+2520*p118+12600*p424-2520*p613-120*p703-12600*p325+840*p604-2520*p505+4200*p406+4200*p334+840*p136-2520*p127)*v+2520*p127)*v)*v+(-2520*p505
  //      +3150*p406-360*p109-2520*p307+1260*p604+45*p802+45*p0010-360*p703+1260*p208+(-720*p118+2520*p217+90*p019+2520*p613-90*p0010+5040*p505-6300*p406-720*p712
  //      -2520*p208-2520*p604+90*p811-5040*p514+720*p109+720*p703-5040*p316-90*p802+6300*p415+5040*p307+(5040*p514-6300*p415+5040*p316-2520*p217-2520*p307+1260*p208
  //      -360*p109-90*p019+1260*p226+1260*p622+720*p712+45*p028+45*p820+45*p0010-2520*p523-360*p721+720*p118-90*p811+3150*p424-2520*p613+45*p802-360*p703-2520*p325
  //      +1260*p604-2520*p505+3150*p406-360*p127)*v)*v+(1260*p505+840*p307-1260*p406-360*p208+10*p901+90*p109-90*p802+360*p703-10*p0010-840*p604+(-10*p019-840*p613
  //      -90*p109-90*p811-1260*p505-360*p217+1260*p406+360*p208+360*p712-840*p307-1260*p415+90*p802+10*p910+90*p118+1260*p514+10*p0010-10*p901-360*p703+840*p316
  //      +840*p604)*v+(-10*p901-10*p109-252*p505+p1000-120*p703+45*p208-120*p307+45*p802+210*p604+210*p406+p0010)*u)*u)*u)*u)*u)*u)*u)*u)*u)*u);
}

void main() {

    // pcaplan
    FragColor = vec4(1,2,3,4);
    return;
    // end
     vec3  rgb;
     vec3  dist;
     float mixVal, kc;
     int  i = 0, idx =0, idxIso = 0;
     vec4 color;


      if ( Picking ) {
       FragColor = encode_id(BasePicking+gl_PrimitiveID);
       return;
      }

     if ( SolOn && Sol == 1 ) {

       //--- Get real position
       vec3 pos = XYZPosition/scenescaling.w + scenescaling.xyz;

        float sol = 0;

       if ( solDeg == 0)
         sol = solp0();
       else  if (  solDeg == 1 )
         sol = solp1(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 2 )
         sol = solp2(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 3 )
         sol = solp3(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 4 )
         sol = solp4(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 5 )
         sol = solp5(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 6 )
         sol = solp6(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 7 )
         sol = solp7(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 8 )
         sol = solp8(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 9 )
         sol = solp9(GBarycentric.x,GBarycentric.y);
       else if (  solDeg == 10 )
         sol = solp10(GBarycentric.x,GBarycentric.y);

       float kc = 0.0;
       if ( sol <= Palette[0] ) {
         sol = Palette[0];
         idx = 1;
       }
       else if ( sol >= Palette[4] ) {
         sol = Palette[4];
         idx = 4;
       }
       else {
         for(i=1; i<5; i++) {
           if ( sol >= Palette[i-1] )
             idx = i;
         }
       }
       for (i=0; i<9; i++) {
          if ( sol >= float(9 - i)/float(9)*Palette[idx - 1] + float(i)/float(9)*Palette[idx] )
            idxIso = i;
       }
       //idxIso = 0;

      kc = (sol - Palette[idx-1]) / (Palette[idx] - Palette[idx-1]);

        rgb = colormap_rgb(idx - 1 + kc, 0, 4);

        //-- isovalues
        float dsol = fwidth(sol);
      //float d2sol = fwidth(dsol);
      //dsol = dsol + 0.5*sqrt(d2sol*min(abs( sol - (float(9 - idxIso)/float(9)*Palette[idx - 1] + float(idxIso)/float(9)*Palette[idx]) ),abs( sol - (float(9 - idxIso - 1)/float(9)*Palette[idx - 1] + float(idxIso + 1)/float(9)*Palette[idx]) )));

        if ( IsoOn ) {
          if ( ( ( abs( sol - (float(9 - idxIso)/float(9)*Palette[idx - 1] + float(idxIso)/float(9)*Palette[idx]) ) < 0.5*IsoWidth*dsol ) )|| (  abs( sol - (float(9 - idxIso - 1)/float(9)*Palette[idx - 1] + float(idxIso + 1)/float(9)*Palette[idx]) ) < 0.5*IsoWidth*dsol) )
            rgb = vec3(0., 0., 0.);
        }

       if( gl_FrontFacing ) {
         if (      UsePhong == 0 )   color = vec4( diffuseModelSol(GPosition,  GNormal , rgb), 1.0 );
         else if ( UsePhong == 1 )   color = vec4( phongModelSol(GPosition,  GNormal , rgb), 1.0 );
         else if ( UsePhong == 2 )   color = vec4( toonModelSol(GNormal, rgb), 1.0 );
         else if ( UsePhong == 3 )   color = vec4( rgb, 1.0 );

       }
       else {
         if (      UsePhong == 0 )   color = vec4( diffuseModelSol(GPosition, -GNormal , rgb  ), 1.0 );
         else if ( UsePhong == 1 )   color = vec4( phongModelSol(GPosition,   -GNormal , rgb  ), 1.0 );
         else if ( UsePhong == 2 )   color = vec4( toonModelSol(-GNormal, rgb), 1.0 );
         else if ( UsePhong == 3 )   color = vec4( rgb, 1.0 );
       }

    }
    else {
      if( gl_FrontFacing ) {
        if (      UsePhong == 0 )   color = vec4( diffuseModel(GPosition,  GNormal   ), 1.0 );
        else if ( UsePhong == 1 )   color = vec4( phongModel(GPosition,  GNormal   ), 1.0 );
        else if ( UsePhong == 2 )   color = vec4( toonModel(GNormal), 1.0 );
        else if ( UsePhong == 3 )   color = vec4( Material.Kd, 1.0 );
        //else if ( UsePhong == 3 )   color = vec4(0.,1.,0.,1.); //Plot of analatycal solution

        if ( CheckInOut    )        color = mix(color,vec4(1.0,0.0,0.0,1.),0.7);
      }
      else {
        if (      UsePhong == 0 )   color = vec4( diffuseModel(GPosition, -GNormal   ), 1.0 );
        else if ( UsePhong == 1 )   color = vec4( phongModel(GPosition,   -GNormal   ), 1.0 );
        else if ( UsePhong == 2 )   color = vec4( toonModel(-GNormal), 1.0 );
        else if ( UsePhong == 3 )   color = vec4( Material.Kd, 1.0 );
        //else if ( UsePhong == 3 )   color = vec4(0.,1.,0.,1.); //Plot of analatycal solution

        if ( CheckInOut    )        color = mix(color,vec4(0.0,0.0,1.0,1.),0.7);
      }
    }


    bool find = false;
    for(int i=0; i<NbrHighLight; i++) {
      if ( gl_PrimitiveID == IdHighlight[i] ) {
        color = mix(color,vec4(1.0,0.0,0.0,1.),0.7);
        color = vec4(0,255,255,1);
        find = true;
      }
    }


    if ( !find && DiscardUnPicked ) discard;

    if ( HighLightAll ) color = mix(color,vec4(1.0,0.0,0.0,1.),0.6);

    if ( !LineOn && !TessOn) { FragColor = color; return; }


   // Find the smallest distance
    float d;
    if ( TessOn )
     dist = GTessDistance;
    else
      dist = GEdgeDistance;
    //  -- compute the norm of the derivative
    vec3 deltas = fwidth(dist);
    dist = vec3( dist.x/deltas.x, dist.y/deltas.y, dist.z/deltas.z  );



     d = min( dist.x, dist.y );
     d = min( d, dist.z );

   if( d < Line.Width - 1 ) {
       mixVal = 1.0;
   } else if( d > Line.Width + 1 ) {
       mixVal = 0.0;
       if ( FaceOff ) discard;
   } else {
       float x = d - (Line.Width - 1);
       mixVal = exp2(-2.0 * (x*x));
   }


   FragColor = mix( color, Line.Color, mixVal );



}
)"
