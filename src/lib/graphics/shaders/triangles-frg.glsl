#version 410


layout( location = 0 ) out vec4 fragColor;

in vec2 v_Parameter;

uniform samplerBuffer solution;
uniform samplerBuffer colormap;

uniform int nb_basis;

// TODO: make these uniforms
const int ncolor = 256;
const float umin = 0;
const float umax = 1;

void
get_color( float u , out vec3 color ) {

    int indx = int(ncolor*(u - umin)/(umax - umin));

    float r0 = texelFetch( colormap , 3*(indx+1) + 0 ).x;
    float g0 = texelFetch( colormap , 3*(indx+1) + 1 ).x;
    float b0 = texelFetch( colormap , 3*(indx+1) + 2 ).x;

    color = vec3(r0,g0,b0);
}

void main() {

    float s = v_Parameter.x;
    float t = v_Parameter.y;

    int idx  = gl_PrimitiveID*nb_basis;
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

    float f = 0.0;
    #if 0
    float phi0 =  s*-3.0-t*3.0+s*t*4.0+(s*s)*2.0+(t*t)*2.0+1.0;
    float phi1 =  -s+(s*s)*2.0;
    float phi2 =  -t+(t*t)*2.0;
    float phi3 =  s*t*4.0;
    float phi4 =  t*(s+t-1.0)*-4.0;
    float phi5 =  -s*(s*4.0+t*4.0-4.0);
    f = f0*phi0 + f1*phi1 + f2*phi2 + f3*phi3 + f4*phi4 + f5*phi5;
    #elif 0

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

    //for (int i = 0; i < 10; i++)
    //    f += u[i]*phi[i];
    f = u[0]*phi[0] + u[1]*phi[1] + u[2]*phi[2] + u[3]*phi[3] + u[4]*phi[4] + u[5]*phi[5] + u[6]*phi[6] + u[7]*phi[7] + u[8]*phi[8] + u[9]*phi[9];

    #elif 0
    float phi0 =  s*(-1.1E1/2.0)-t*(1.1E1/2.0)+s*t*1.8E1-s*(t*t)*(2.7E1/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*9.0-
                (s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*t*t)*(9.0/2.0)+1.0;

    float phi1 =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);

    float phi2 =  t-(t*t)*(9.0/2.0)+(t*t*t)*(9.0/2.0);

    float phi3 =  s*t*(-9.0/2.0)+(s*s)*t*(2.7E1/2.0);

    float phi4 =  s*t*(-9.0/2.0)+s*(t*t)*(2.7E1/2.0);

    float phi5 =  t*(-9.0/2.0)+s*t*(9.0/2.0)-s*(t*t)*(2.7E1/2.0)+(t*t)*1.8E1-(t*t*t)*(2.7E1/2.0);

    float phi6 =  t*9.0-s*t*(4.5E1/2.0)+s*(t*t)*2.7E1+(s*s)*t*(2.7E1/2.0)-(t*t)*(4.5E1/2.0)+(t*t*t)*(2.7E1/2.0);

    float phi7 =  s*9.0-s*t*(4.5E1/2.0)+s*(t*t)*(2.7E1/2.0)+(s*s)*t*2.7E1-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0);

    float phi8 =  s*(-9.0/2.0)+s*t*(9.0/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*1.8E1-(s*s*s)*(2.7E1/2.0);

    float phi9 =  s*t*2.7E1-s*(t*t)*2.7E1-(s*s)*t*2.7E1;

    f = phi0*f0 + phi1*f1 + phi2*f2 + phi3*f3 + phi4*f4 + phi5*f5 + phi6*f6 + phi7*f7 + phi8*f8 + phi9*f9;
    #elif 0
    f =  f0*(s*(-5.5)-t*(5.5)+s*t*18-s*(t*t)*(13.5)-(s*s)*t*(13.5)+(s*s)*9.0-
                (s*s*s)*(4.5)+(t*t)*9.0-(t*t*t)*(4.5)+1.0)

                + f1*(s-(s*s)*(4.5)+(s*s*s)*(4.5))

                + f2*(t-(t*t)*(4.5)+(t*t*t)*(4.5))

                + f3*(s*t*(-4.5)+(s*s)*t*(13.5))

                + f4*(s*t*(-9.0/2.0)+s*(t*t)*(13.5))

                + f5*(t*(-9.0/2.0)+s*t*(4.5)-s*(t*t)*(13.5)+(t*t)*18-(t*t*t)*(13.5))

                + f6*(t*9.0-s*t*(22.5)+s*(t*t)*27+(s*s)*t*(13.5)-(t*t)*(22.5)+(t*t*t)*(13.5))

                + f7*(s*9.0-s*t*(22.5)+s*(t*t)*(13.5)+(s*s)*t*27-(s*s)*(22.5)+(s*s*s)*(13.5))

                + f8*(s*(-4.5)+s*t*(4.5)-(s*s)*t*(13.5)+(s*s)*18-(s*s*s)*(13.5))

                + f9*(s*t*27-s*(t*t)*27-(s*s)*t*27);
    #endif

    vec3 color;
    get_color(s+t,color);
    //color = vec3(0.2,0.8,0.2);
    fragColor = vec4(color,1.0);
}
