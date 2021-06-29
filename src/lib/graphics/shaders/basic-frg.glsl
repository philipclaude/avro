#version 400


layout( location = 0 ) out vec4 fragColor;

in vec3 v_Position;
noperspective in vec2 v_Parameter;

uniform samplerBuffer solution;
//uniform samplerBuffer colormap;

int
get_color_idx( float u ) {
    float umin = 0.3;
    float umax = 0.7;

    float frac = 255.0*(u - umin)/(umax - umin);
    return int(frac);
}

void main() {

    int idx = 0;//gl_PrimitiveID*6;

    float u = v_Parameter.x;
    float v = v_Parameter.y;

    float f0 = texelFetch( solution , idx + 0 ).x;
    float f1 = texelFetch( solution , idx + 1 ).x;
    float f2 = texelFetch( solution , idx + 2 ).x;
    //float f3 = texelFetch( solution , idx + 3 ).x;
    //float f4 = texelFetch( solution , idx + 4 ).x;
    //float f5 = texelFetch( solution , idx + 5 ).x;

    float f = f0 * u + f1 * v + f2 * (1 - u - v);
    int color_idx = get_color_idx(f);

    /*
    float r = texelFetch( colormap , 3*color_idx + 0 ).x;
    float g = texelFetch( colormap , 3*color_idx + 1 ).x;
    float b = texelFetch( colormap , 3*color_idx + 2 ).x;
    */

    fragColor = vec4(f,f,f,1.0);
}
