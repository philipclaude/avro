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


uniform MaterialInfo Material;

uniform struct LineInfo {
  float Width;
  vec4 Color;
} Line;

//in vec3 GPosition;

uniform int  NbrHighLight;
uniform int  IdHighlight[10];
uniform bool HighLightAll;


/* display option */
uniform int  BasePicking;
uniform bool DiscardUnPicked;
uniform bool LineOn;
uniform bool CheckInOut;
uniform bool Picking;
uniform int  UsePhong;
uniform int  TypePrim;
uniform bool FaceOff;



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
    //for debug
    //if ( id == 0 ) return vec4( 1, 1, 1, 1.0);
    //if ( id == 1 ) return vec4( 1, 0, 0, 1.0);
    //if ( id == 2 ) return vec4( 0, 1, 0, 1.0);
    //if ( id == 3 ) return vec4( 0, 0, 1, 1.0);
    //if ( id == 4 ) return vec4( 1, 0, 1, 1.0);
    //if ( id == 5 ) return vec4( 1, 0, 1, 1.0);
    //if ( id == 6 ) return vec4( 1, 0, 1, 1.0);
    //if ( id == 4 ) return vec4( 1, 1, 1, 1.0);
    //if ( id == 5 ) return vec4( 1, 1, 1, 1.0);
    //if ( id == 6 ) return vec4( 1, 1, 1, 1.0);


	// convert to floats. only divide by 255, because range is 0-255
	return color;
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

void main() {

     vec4 color;

     if ( Picking ) {
       FragColor = encode_id(BasePicking+gl_PrimitiveID);
       return;
     }

    color = vec4( Material.Kd, 1.0 );

    //FragColor = color;
    //FragColor = vec4( vec3(color) ,1);

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

    FragColor = color;

    return;
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

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

//#define resolution vec2(50.0, 50.0)
//#define Thickness 3
//
//float drawLine(vec2 p1, vec2 p2) {
//  vec2 uv = gl_FragCoord.xy / resolution.xy;
//
//  float a = abs(distance(p1, uv));
//  float b = abs(distance(p2, uv));
//  float c = abs(distance(p1, p2));
//
//  if ( a >= c || b >=  c ) return 0.0;
//
//  float p = (a + b + c) * 0.5;
//
//  // median to (p1, p2) vector
//  float h = 2 / c * sqrt( p * ( p - a) * ( p - b) * ( p - c));
//
//  return mix(1.0, 0.0, smoothstep(0.5 * Thickness, 1.5 * Thickness, h));
//}

//void main()
//{
//  FragColor = vec4(
//      max(
//        max(
//          drawLine(vec2(0.1, 0.1), vec2(0.1, 0.9)),
//          drawLine(vec2(0.1, 0.9), vec2(0.7, 0.5))),
//          drawLine(vec2(0.1, 0.1), vec2(0.7, 0.5))));
//}

)"
