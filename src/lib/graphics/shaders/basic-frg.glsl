#version 330


layout( location = 0 ) out vec4 fragColor;

in vec3 v_Position;
in vec2 v_Parameter;

uniform samplerBuffer solution;
uniform samplerBuffer colormap;

vec3
get_color( float u ) {

    int ncolor = 256;

    float umin = -4.0;
    float umax =  4.0;

    float frac = ncolor*(u - umin)/(umax - umin);
    int   indx = int(frac);
    frac -= float(indx);

    float r0 = texelFetch( colormap , 3*(indx+1) + 0 ).x;
    float g0 = texelFetch( colormap , 3*(indx+1) + 1 ).x;
    float b0 = texelFetch( colormap , 3*(indx+1) + 2 ).x;

    float r1 = texelFetch( colormap , 3*(indx  ) + 0 ).x;
    float g1 = texelFetch( colormap , 3*(indx  ) + 1 ).x;
    float b1 = texelFetch( colormap , 3*(indx  ) + 2 ).x;

    vec3 c0 = vec3(r0,g0,b0);
    vec3 c1 = vec3(r1,g1,b1);

    return frac*c0 + (1 - frac)*c1;
}

void main() {

    int idx = 0;//gl_PrimitiveID*6;

    float u = v_Parameter.x;
    float v = v_Parameter.y;

    float f0 = texelFetch( solution , idx + 0 ).x;
    float f1 = texelFetch( solution , idx + 1 ).x;
    float f2 = texelFetch( solution , idx + 2 ).x;
    float f3 = texelFetch( solution , idx + 3 ).x;
    float f4 = texelFetch( solution , idx + 4 ).x;
    float f5 = texelFetch( solution , idx + 5 ).x;

    float s = u, t = v;
    float phi0 =  s*-3.0-t*3.0+s*t*4.0+(s*s)*2.0+(t*t)*2.0+1.0;
    float phi1 =  -s+(s*s)*2.0;
    float phi2 =  -t+(t*t)*2.0;
    float phi3 =  s*t*4.0;
    float phi4 =  t*(s+t-1.0)*-4.0;
    float phi5 =  -s*(s*4.0+t*4.0-4.0);

    float f = f0 * (1 - u - v) + f1 * u + f2 * v;
    f = f0*phi0 + f1*phi1 + f2*phi2 + f3*phi3 + f4*phi4 + f5*phi5;
    vec3 color = get_color(f);

    fragColor = vec4(color,1.0);
}
