#include "common/tools.h"

#include "library/eps.h"

#include <cstdlib>

namespace luna
{

namespace library
{

epsFile::epsFile() :
  fid_(NULL),
  viewport_(4,0)
{}

void
epsFile::write( const std::string& filename )
{
  fid_ = fopen(filename.c_str(),"w");

  fprintf(fid_,"%%!PS-Adobe-3.0 EPSF-3.0\n");

  fprintf(fid_,"%%%%Title: some plot\n"
              "%%%%Creator: luna 1.0.0b, (c) Philip Claude Caplan\n"
              "%%%%For: some reason\n"
              "%%%%CreationDate: October 27 2019\n"
              "%%%%LanguageLevel: 3\n"
              "%%%%DocumentData: Clean7Bit\n"
              "%%%%Pages: 1\n");

   fprintf(fid_,"%%%%BoundingBox: %d %d %d %d\n"
              "%%%%EndComments\n",
              viewport_[0],
              viewport_[1],
              viewport_[2],
              viewport_[3] );

    float threshold[3] = {0.064F,0.034F,0.100F};

    fprintf(fid_,"%%%%BeginProlog\n"
                "/gl2psdict 64 dict def gl2psdict begin\n"
                "/tryPS3shading %s def %% set to false to force subdivision\n"
                "/rThreshold %g def %% red component subdivision threshold\n"
                "/gThreshold %g def %% green component subdivision threshold\n"
                "/bThreshold %g def %% blue component subdivision threshold\n",
                //(gl2ps->options & GL2PS_NO_PS3_SHADING) ? "false" : "true",
                "true",
                threshold[0], threshold[1], threshold[2]);

    fprintf(fid_,"/BD { bind def } bind def\n"
                "/C  { setrgbcolor } BD\n"
                "/G  { 0.082 mul exch 0.6094 mul add exch 0.3086 mul add neg 1.0 add setgray } BD\n"
                "/W  { setlinewidth } BD\n"
                "/LC  { setlinecap } BD\n"
                "/LJ  { setlinejoin } BD\n");

    fprintf(fid_,"/FC { findfont exch /SH exch def SH scalefont setfont } BD\n"
                "/SW { dup stringwidth pop } BD\n"
                "/S  { FC moveto show } BD\n"
                "/SBC{ FC moveto SW -2 div 0 rmoveto show } BD\n"
                "/SBR{ FC moveto SW neg 0 rmoveto show } BD\n"
                "/SCL{ FC moveto 0 SH -2 div rmoveto show } BD\n"
                "/SCC{ FC moveto SW -2 div SH -2 div rmoveto show } BD\n"
                "/SCR{ FC moveto SW neg SH -2 div rmoveto show } BD\n"
                "/STL{ FC moveto 0 SH neg rmoveto show } BD\n"
                "/STC{ FC moveto SW -2 div SH neg rmoveto show } BD\n"
                "/STR{ FC moveto SW neg SH neg rmoveto show } BD\n");

    /* rotated text routines: same nameanem with R appended */

    fprintf(fid_,"/FCT { FC translate 0 0 } BD\n"
                "/SR  { gsave FCT moveto rotate show grestore } BD\n"
                "/SBCR{ gsave FCT moveto rotate SW -2 div 0 rmoveto show grestore } BD\n"
                "/SBRR{ gsave FCT moveto rotate SW neg 0 rmoveto show grestore } BD\n"
                "/SCLR{ gsave FCT moveto rotate 0 SH -2 div rmoveto show grestore} BD\n");
    fprintf(fid_,"/SCCR{ gsave FCT moveto rotate SW -2 div SH -2 div rmoveto show grestore} BD\n"
                "/SCRR{ gsave FCT moveto rotate SW neg SH -2 div rmoveto show grestore} BD\n"
                "/STLR{ gsave FCT moveto rotate 0 SH neg rmoveto show grestore } BD\n"
                "/STCR{ gsave FCT moveto rotate SW -2 div SH neg rmoveto show grestore } BD\n"
                "/STRR{ gsave FCT moveto rotate SW neg SH neg rmoveto show grestore } BD\n");

    fprintf(fid_,"/P  { newpath 0.0 360.0 arc closepath fill } BD\n"
                "/LS { newpath moveto } BD\n"
                "/L  { lineto } BD\n"
                "/LE { lineto stroke } BD\n"
                "/T  { newpath moveto lineto lineto closepath fill } BD\n");

    /* Smooth-shaded triangle with PostScript level 3 shfill operator:
          x3 y3 r3 g3 b3 x2 y2 r2 g2 b2 x1 y1 r1 g1 b1 STshfill */

    fprintf(fid_,"/STshfill {\n"
                "      /b1 exch def /g1 exch def /r1 exch def /y1 exch def /x1 exch def\n"
                "      /b2 exch def /g2 exch def /r2 exch def /y2 exch def /x2 exch def\n"
                "      /b3 exch def /g3 exch def /r3 exch def /y3 exch def /x3 exch def\n"
                "      gsave << /ShadingType 4 /ColorSpace [/DeviceRGB]\n"
                "      /DataSource [ 0 x1 y1 r1 g1 b1 0 x2 y2 r2 g2 b2 0 x3 y3 r3 g3 b3 ] >>\n"
                "      shfill grestore } BD\n");

    /* Flat-shaded triangle with middle color:
          x3 y3 r3 g3 b3 x2 y2 r2 g2 b2 x1 y1 r1 g1 b1 Tm */

    fprintf(fid_,/* stack : x3 y3 r3 g3 b3 x2 y2 r2 g2 b2 x1 y1 r1 g1 b1 */
                "/Tm { 3 -1 roll 8 -1 roll 13 -1 roll add add 3 div\n" /* r = (r1+r2+r3)/3 */
                /* stack : x3 y3 g3 b3 x2 y2 g2 b2 x1 y1 g1 b1 r */
                "      3 -1 roll 7 -1 roll 11 -1 roll add add 3 div\n" /* g = (g1+g2+g3)/3 */
                /* stack : x3 y3 b3 x2 y2 b2 x1 y1 b1 r g b */
                "      3 -1 roll 6 -1 roll 9 -1 roll add add 3 div" /* b = (b1+b2+b3)/3 */
                /* stack : x3 y3 x2 y2 x1 y1 r g b */
                " C T } BD\n");

    /* Split triangle in four sub-triangles (at sides middle points) and call the
       STnoshfill procedure on each, interpolating the colors in RGB space:
          x3 y3 r3 g3 b3 x2 y2 r2 g2 b2 x1 y1 r1 g1 b1 STsplit
       (in procedure comments key: (Vi) = xi yi ri gi bi) */

    fprintf(fid_,"/STsplit {\n"
                "      4 index 15 index add 0.5 mul\n" /* x13 = (x1+x3)/2 */
                "      4 index 15 index add 0.5 mul\n" /* y13 = (y1+y3)/2 */
                "      4 index 15 index add 0.5 mul\n" /* r13 = (r1+r3)/2 */
                "      4 index 15 index add 0.5 mul\n" /* g13 = (g1+g3)/2 */
                "      4 index 15 index add 0.5 mul\n" /* b13 = (b1+b3)/2 */
                "      5 copy 5 copy 25 15 roll\n");

    /* at this point, stack = (V3) (V13) (V13) (V13) (V2) (V1) */

    fprintf(fid_,"      9 index 30 index add 0.5 mul\n" /* x23 = (x2+x3)/2 */
                "      9 index 30 index add 0.5 mul\n" /* y23 = (y2+y3)/2 */
                "      9 index 30 index add 0.5 mul\n" /* r23 = (r2+r3)/2 */
                "      9 index 30 index add 0.5 mul\n" /* g23 = (g2+g3)/2 */
                "      9 index 30 index add 0.5 mul\n" /* b23 = (b2+b3)/2 */
                "      5 copy 5 copy 35 5 roll 25 5 roll 15 5 roll\n");

    /* stack = (V3) (V13) (V23) (V13) (V23) (V13) (V23) (V2) (V1) */

    fprintf(fid_,"      4 index 10 index add 0.5 mul\n" /* x12 = (x1+x2)/2 */
                "      4 index 10 index add 0.5 mul\n" /* y12 = (y1+y2)/2 */
                "      4 index 10 index add 0.5 mul\n" /* r12 = (r1+r2)/2 */
                "      4 index 10 index add 0.5 mul\n" /* g12 = (g1+g2)/2 */
                "      4 index 10 index add 0.5 mul\n" /* b12 = (b1+b2)/2 */
                "      5 copy 5 copy 40 5 roll 25 5 roll 15 5 roll 25 5 roll\n");

    /* stack = (V3) (V13) (V23) (V13) (V12) (V23) (V13) (V1) (V12) (V23) (V12) (V2) */

    fprintf(fid_,"      STnoshfill STnoshfill STnoshfill STnoshfill } BD\n");

    /* Gouraud shaded triangle using recursive subdivision until the difference
       between corner colors does not exceed the thresholds:
          x3 y3 r3 g3 b3 x2 y2 r2 g2 b2 x1 y1 r1 g1 b1 STnoshfill  */

    fprintf(fid_,"/STnoshfill {\n"
                "      2 index 8 index sub abs rThreshold gt\n" /* |r1-r2|>rth */
                "      { STsplit }\n"
                "      { 1 index 7 index sub abs gThreshold gt\n" /* |g1-g2|>gth */
                "        { STsplit }\n"
                "        { dup 6 index sub abs bThreshold gt\n" /* |b1-b2|>bth */
                "          { STsplit }\n"
                "          { 2 index 13 index sub abs rThreshold gt\n" /* |r1-r3|>rht */
                "            { STsplit }\n"
                "            { 1 index 12 index sub abs gThreshold gt\n" /* |g1-g3|>gth */
                "              { STsplit }\n"
                "              { dup 11 index sub abs bThreshold gt\n" /* |b1-b3|>bth */
                "                { STsplit }\n"
                "                { 7 index 13 index sub abs rThreshold gt\n"); /* |r2-r3|>rht */
    fprintf(fid_,"                  { STsplit }\n"
                "                  { 6 index 12 index sub abs gThreshold gt\n" /* |g2-g3|>gth */
                "                    { STsplit }\n"
                "                    { 5 index 11 index sub abs bThreshold gt\n" /* |b2-b3|>bth */
                "                      { STsplit }\n"
                "                      { Tm }\n" /* all colors sufficiently similar */
                "                      ifelse }\n"
                "                    ifelse }\n"
                "                  ifelse }\n"
                "                ifelse }\n"
                "              ifelse }\n"
                "            ifelse }\n"
                "          ifelse }\n"
                "        ifelse }\n"
                "      ifelse } BD\n");

    fprintf(fid_,"tryPS3shading\n"
                "{ /shfill where\n"
                "  { /ST { STshfill } BD }\n"
                "  { /ST { STnoshfill } BD }\n"
                "  ifelse }\n"
                "{ /ST { STnoshfill } BD }\n"
                "ifelse\n");

    fprintf(fid_,"end\n"
                "%%%%EndProlog\n"
                "%%%%BeginSetup\n"
                "/DeviceRGB setcolorspace\n"
                "gl2psdict begin\n"
                "%%%%EndSetup\n"
                "%%%%Page: 1 1\n"
                "%%%%BeginPageSetup\n");

    fprintf(fid_,"%%%%EndPageSetup\n"
                "mark\n"
                "gsave\n"
                "%1.1f 1.0 scale\n",1024./640.);

  int x = viewport_[0], y = viewport_[1], w = viewport_[2], h = viewport_[3];

  fprintf(fid_,"newpath %d %d moveto %d %d lineto %d %d lineto %d %d lineto\n"
            "closepath clip\n",
            x, y, x+w, y, x+w, y+h, x, y+h);



  print_triangles();

  fprintf(fid_,"grestore\n"
            "showpage\n"
            "cleartomark\n"
            "%%%%PageTrailer\n"
            "%%%%Trailer\n"
            "end\n"
            "%%%%EOF\n");


  fclose(fid_);
}

void
epsFile::set_viewport( int* viewport )
{
  viewport_[0] = viewport[0];
  viewport_[1] = viewport[1];
  viewport_[2] = viewport[2];
  viewport_[3] = viewport[3];
}

void
epsFile::add_triangles( const std::vector<real_t>& t , const std::vector<real_t>& c )
{
  for (index_t k=0;k<t.size();k++)
    triangles_.push_back(100*t[k]);
  for (index_t k=0;k<c.size();k++)
    triangle_colors_.push_back(c[k]);
}

void
get_constant_color( const real_t* data , float* c )
{
  c[0] = ( data[0] + data[3] + data[6] )/3.;
  c[1] = ( data[1] + data[4] + data[7] )/3.;
  c[2] = ( data[2] + data[5] + data[8] )/3.;
}

void
epsFile::print_triangles() const
{
  int w = 640/2,h = 320;
  real_t xnor,ynor,znor;
  std::vector<real_t> zbuffer;
  std::vector<index_t> visible;
  std::vector<real_t> distance( triangles_.size()/9 , -1 );
  std::vector<index_t> idx;
  index_t count = 0;
  for (index_t k=0;k<triangles_.size()/9;k++)
  {
    // compute the triangle normal
    real_t x1 = triangles_[9*k  ];
    real_t y1 = triangles_[9*k+1];
    real_t z1 = triangles_[9*k+2];
    real_t x2 = triangles_[9*k+3];
    real_t y2 = triangles_[9*k+4];
    real_t z2 = triangles_[9*k+5];
    real_t x3 = triangles_[9*k+6];
    real_t y3 = triangles_[9*k+7];
    real_t z3 = triangles_[9*k+8];

    xnor = (y3-y1)*(z2-z1)-(y2-y1)*(z3-z1);
    ynor = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1);
    znor = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);

    if (znor > 0.0) continue;

    idx.push_back(count++);
    visible.push_back( k );

    real_t z0 = (z1+z2+z3)/3.0;
    zbuffer.push_back( z0 );
    distance[k] = z0;
  }

  std::sort( idx.begin() , idx.end() , SortBy<real_t>(zbuffer) );
  std::reverse( idx.begin() , idx.end() );
  for (index_t i=0;i<idx.size();i++)
  {
    index_t k = visible[idx[i]];

    #if 1
    float c[3];
    get_constant_color( &triangle_colors_[9*k] , c );
    fprintf(fid_,"%g %g %g C\n",c[0],c[1],c[2]);
    fprintf(fid_,"%g %g %g %g %g %g T\n",
        triangles_[9*k]  +w,triangles_[9*k+1]+h,
        triangles_[9*k+3]+w,triangles_[9*k+4]+h,
        triangles_[9*k+6]+w,triangles_[9*k+7]+h);
    #else
    fprintf(fid_,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ST\n",
        triangles_[9*k]+w,triangles_[9*k+1]+h,
        triangle_colors_[9*k],triangle_colors_[9*k+1],triangle_colors_[9*k+2],

        triangles_[9*k+3]+w,triangles_[9*k+4]+h,
        triangle_colors_[9*k+3],triangle_colors_[9*k+4],triangle_colors_[9*k+5],

        triangles_[9*k+6]+w,triangles_[9*k+7]+h,
        triangle_colors_[9*k+6],triangle_colors_[9*k+7],triangle_colors_[9*k+8] );
    #endif
  }

}

} // library

} // luna
