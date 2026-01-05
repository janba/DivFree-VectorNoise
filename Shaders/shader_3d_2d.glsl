/*

  Accompanying code, Improving Curl Noise, SIGGRAPH Asia 2025
  ---

  This demonstrates curl noise within a volume, progressively
  constrained to follow the enclosing surface.
  The displacement vector is obtained as the cross product of
  a solid noise gradient and the gradient of another noise that
  is progressively interpolated to become the enclosing surface
  normal.

*/

// -------------------------------------------------------
// -------------------------------------------------------
// Support code for noise and sphere rendering, skip over
// 'end of support code' to reach our method part.
// (We use code from Inigo Quilez).
// -------------------------------------------------------
// -------------------------------------------------------

// The MIT License
// Copyright © 2017 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// https://www.youtube.com/c/InigoQuilez
// https://iquilezles.org/
//
// Computes the analytic derivatives of a 3D Gradient Noise. This can be used for example to compute normals to a
// 3d rocks based on Gradient Noise without approximating the gradient by having to take central differences.
//
// More info here: https://iquilezles.org/articles/gradientnoise
//
vec3 hash( ivec3 p )     // this hash is not production ready, please
{                        // replace this by something better
	ivec3 n = ivec3( p.x*127 + p.y*311 + p.z*74,
                     p.x*269 + p.y*183 + p.z*246,
                     p.x*113 + p.y*271 + p.z*124);

	// 1D hash by Hugo Elias
	n = (n << 13) ^ n;
    n = n * (n * n * 15731 + 789221) + 1376312589;
    return -1.0+2.0*vec3( n & ivec3(0x0fffffff))/float(0x0fffffff);
}

// return value noise (in x) and its derivatives (in yzw)
vec4 noised( in vec3 x )
{
    // grid
    ivec3 i = ivec3(floor(x));
    vec3 f = fract(x);
    // quintic interpolant
    vec3 u = f*f*f*(f*(f*6.0-15.0)+10.0);
    vec3 du = 30.0*f*f*(f*(f-2.0)+1.0);
    // gradients
    vec3 ga = hash( i+ivec3(0,0,0) );
    vec3 gb = hash( i+ivec3(1,0,0) );
    vec3 gc = hash( i+ivec3(0,1,0) );
    vec3 gd = hash( i+ivec3(1,1,0) );
    vec3 ge = hash( i+ivec3(0,0,1) );
	vec3 gf = hash( i+ivec3(1,0,1) );
    vec3 gg = hash( i+ivec3(0,1,1) );
    vec3 gh = hash( i+ivec3(1,1,1) );
    // projections
    float va = dot( ga, f-vec3(0.0,0.0,0.0) );
    float vb = dot( gb, f-vec3(1.0,0.0,0.0) );
    float vc = dot( gc, f-vec3(0.0,1.0,0.0) );
    float vd = dot( gd, f-vec3(1.0,1.0,0.0) );
    float ve = dot( ge, f-vec3(0.0,0.0,1.0) );
    float vf = dot( gf, f-vec3(1.0,0.0,1.0) );
    float vg = dot( gg, f-vec3(0.0,1.0,1.0) );
    float vh = dot( gh, f-vec3(1.0,1.0,1.0) );
    // interpolations
    return vec4( va + u.x*(vb-va) + u.y*(vc-va) + u.z*(ve-va) + u.x*u.y*(va-vb-vc+vd) + u.y*u.z*(va-vc-ve+vg) + u.z*u.x*(va-vb-ve+vf) + (-va+vb+vc-vd+ve-vf-vg+vh)*u.x*u.y*u.z,    // value
                 ga + u.x*(gb-ga) + u.y*(gc-ga) + u.z*(ge-ga) + u.x*u.y*(ga-gb-gc+gd) + u.y*u.z*(ga-gc-ge+gg) + u.z*u.x*(ga-gb-ge+gf) + (-ga+gb+gc-gd+ge-gf-gg+gh)*u.x*u.y*u.z +   // derivatives
                 du * (vec3(vb,vc,ve) - va + u.yzx*vec3(va-vb-vc+vd,va-vc-ve+vg,va-vb-ve+vf) + u.zxy*vec3(va-vb-ve+vf,va-vb-vc+vd,va-vc-ve+vg) + u.yzx*u.zxy*(-va+vb+vc-vd+ve-vf-vg+vh) ));
}

// https://iquilezles.org/articles/intersectors/
float iSphere( in vec3 ro, in vec3 rd, float sph )
{
    float ra = sph;
    vec3  oc = ro;
    float b  = dot( oc, rd );
    float c  = dot( oc, oc ) - ra*ra;
    float h  = b*b - c;
    if( h < 0.0 ) return (-1.0); // no intersection
    h        = sqrt( h );
    return -b-h;
}

// df(x)/dx
vec3 nSphere( in vec3 pos, float sph )
{
    return normalize(pos);
}

float dSphere( vec3 p, float sph )
{
  float ra = sph;
  return length(p)-ra;
}

// -------------------------------------------------------
// -------------------------------------------------------
// end of support code
// -------------------------------------------------------
// -------------------------------------------------------

float pattern(vec3 p)
{
  p = normalize(p);
  float theta = acos(p.z);
  float phi   = atan(p.y, p.x);
  vec2 f = fract(vec2(theta,phi)*3.0);
  vec2 d = 3.0*abs(0.5-f);
  return min(d.x,d.y);
}

const float noise_scale = 4.0;
const float crust = 0.5;

vec3 nVol(vec3 pos,float sph)
{
  float d = dSphere( pos, sph);
  vec3 nS = nSphere( pos, sph );               // on-surface curl noise
  vec3 nN = noised( pos.zyx*noise_scale ).yzw; // 3D noise gradient
  //
  // return nN; // test 3D noise only
  // return nS; // test surface constrained only
  //
  float i = min(1.0, -2.0*d/crust); // interpolation based on sphere SDF
  return nS * (1.0-i) + nN * i;
}

vec3 perturb3(vec3 pos, float sph)
{
    // ----------
    // amount of displacement
    float d = 0.35*sin(iTime);
    // number of steps
    const int N = 16;
    // integrate displacement
    float dd = d / float(N);
    vec3  cn;
    vec4  n4  = noised( pos*noise_scale );
    vec3  nor = nVol( pos, sph );
    float tau = n4.x; // iso-contour value
    for (int i=0;i<N;++i) {
      cn  = cross(nor,n4.yzw);
      pos = pos + cn * dd;
      n4  = noised( pos*noise_scale );
      nor = nVol( pos, sph );
      // project on noise iso-contour
      vec3 g = n4.yzw*noise_scale;
      pos = pos - g * (n4.x-tau) / (dot(g,g)+1e-30);
    }
    return pos;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // camera movement
	float an = iTime*0.5;
	vec3 ro  = 1.1*vec3( 0.9*cos(an), 0.7*sin(an), 1.5 );
    vec3 ta  = vec3( 0.0, 0.0, 0.0 );
    // camera matrix
    vec3 ww  = normalize( ta - ro );
    vec3 uu  = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
    vec3 vv  = normalize( cross(uu,ww));
    vec3 tot = vec3(0.0);
    // pixel coordinates
    vec2 p   = (-iResolution.xy + 2.0*fragCoord)/iResolution.y;
    // create view ray
    vec3 rd  = normalize( p.x*uu + p.y*vv + 1.5*ww );
    // raytrace sphere
    float sph = 1.0;
    float t  = iSphere( ro, rd, sph );
    // shading/lighting
    vec3 col = vec3(0.0);
    if( t>0.0 ) {
        // ray-march
        float density = 0.0;
        bool first = true;
        float ds = 0.0;
        for (float tt = t ; tt < t + 1.0 && density < 1.0 && -ds < crust; tt += 0.0005) {
            vec3  pos = ro + tt*rd;
            ds  = dSphere( pos, sph );
            if (ds > 1e-6) break;
            vec3  nor = nSphere( pos, sph );
            pos       = perturb3( pos,sph );
            float ptn = pattern(pos);
            ptn = ptn < 0.1 ? 1.0 : 0.0;
            if (ds < -1e-6) {
            if (first) {
                density  = 1.0 * ptn;
                // col      = density * vec3(0.4,0.4,0.7) * max(0.0,nor.x);
                col      = density * vec3(0.4,0.7,1.2)
                                   * max(0.0,0.7*nor.x+0.3*nor.y);
                first    = false;
            } else {
                // break;
                float dn = ptn * (-ds < crust ? 1.0 : 0.0);
                col     += vec3(0.25 + max(0.0,nor.x))
                         * abs(sin(-ds*120.0))
                         * dn * (1.0-density)
                         * (1.0+ds*2.0);
                density  = density + (1.0-density) * dn;
            }
            }
        }
        // done
        col += col + (1.0-density) * vec3(1,1,1);
    } else {
        // no hit, background color
        col = vec3(1.0);
    }

    col = sqrt( col );

    tot += col;

	fragColor = vec4( tot, 1.0 );
}
