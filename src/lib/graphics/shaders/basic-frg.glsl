#version 330


layout( location = 0 ) out vec4 fragColor;

in vec3 v_Position;
in vec2 v_Parameter;

uniform samplerBuffer solution;
uniform samplerBuffer colormap;

uniform int nb_basis;

vec3
get_color( float u ) {

    int ncolor = 256;

    float umin = 0;
    float umax =  1;

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

    int idx = gl_PrimitiveID*nb_basis;

    float s = v_Parameter.x;
    float t = v_Parameter.y;

    float f0 = texelFetch( solution , idx + 0 ).x;
    float f1 = texelFetch( solution , idx + 1 ).x;
    float f2 = texelFetch( solution , idx + 2 ).x;
    float f3 = texelFetch( solution , idx + 3 ).x;
    float f4 = texelFetch( solution , idx + 4 ).x;
    float f5 = texelFetch( solution , idx + 5 ).x;

    float f6 = texelFetch( solution , idx + 6 ).x;
    float f7 = texelFetch( solution , idx + 7 ).x;
    float f8 = texelFetch( solution , idx + 8 ).x;
    float f9 = texelFetch( solution , idx + 9 ).x;

    #if 0
    float phi0 =  s*-3.0-t*3.0+s*t*4.0+(s*s)*2.0+(t*t)*2.0+1.0;
    float phi1 =  -s+(s*s)*2.0;
    float phi2 =  -t+(t*t)*2.0;
    float phi3 =  s*t*4.0;
    float phi4 =  t*(s+t-1.0)*-4.0;
    float phi5 =  -s*(s*4.0+t*4.0-4.0);
    float f = f0*phi0 + f1*phi1 + f2*phi2 + f3*phi3 + f4*phi4 + f5*phi5;
    #else

    float u[10];
    u[0] = f0; u[1] = f1; u[2] = f2; u[3] = f3; u[4] = f4; u[5] = f5; u[6] = f6; u[7] = f7; u[8] = f8; u[9] = f9;
    float phi[10];
    phi[0] =  s*(-1.1E1/2.0)-t*(1.1E1/2.0)+s*t*1.8E1-s*(t*t)*(2.7E1/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*9.0-
                (s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*t*t)*(9.0/2.0)+1.0;

    phi[1] =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);

    phi[2] =  t-(t*t)*(9.0/2.0)+(t*t*t)*(9.0/2.0);

    phi[3] =  s*t*(-9.0/2.0)+(s*s)*t*(2.7E1/2.0);

    phi[4] =  s*t*(-9.0/2.0)+s*(t*t)*(2.7E1/2.0);

    phi[5] =  t*(-9.0/2.0)+s*t*(9.0/2.0)-s*(t*t)*(2.7E1/2.0)+(t*t)*1.8E1-(t*t*t)*(2.7E1/2.0);

    phi[6] =  t*9.0-s*t*(4.5E1/2.0)+s*(t*t)*2.7E1+(s*s)*t*(2.7E1/2.0)-(t*t)*(4.5E1/2.0)+(t*t*t)*(2.7E1/2.0);

    phi[7] =  s*9.0-s*t*(4.5E1/2.0)+s*(t*t)*(2.7E1/2.0)+(s*s)*t*2.7E1-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0);

    phi[8] =  s*(-9.0/2.0)+s*t*(9.0/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*1.8E1-(s*s*s)*(2.7E1/2.0);

    phi[9] =  s*t*2.7E1-s*(t*t)*2.7E1-(s*s)*t*2.7E1;

    float f = 0.0;
    for (int i = 0; i < 10; i++)
        f += u[i]*phi[i];
    #endif

    vec3 color = vec3(0,0,0);
    //vec3 color = get_color(f);

    fragColor = vec4(color,1.0);
}
