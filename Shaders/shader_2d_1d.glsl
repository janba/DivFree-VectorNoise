/*

  Accompanying code, Improving Curl Noise, SIGGRAPH Asia 2025
  ---

  This demonstrates curl noise within an area progressively
  constrained to follow the boundary curve.
  The displacement vector is obtained as the cross product of a
  solid noise gradient and a vector interpolated from (0,0,1)
  inside to a circle radial vector outside.

*/

// -------------------------------------------------------
// -------------------------------------------------------
// Support code for noise, skip over
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

// -------------------------------------------------------
// -------------------------------------------------------
// end of support code
// -------------------------------------------------------
// -------------------------------------------------------

const float R = 1.0;

float pattern(vec2 p)
{
  float a = atan(p.y,p.x);
  float r = length(p);
  return fract(a*5.0) > 0.5 ? (0.5+fract(r*2.0)) : 0.0;
}

const float clamp_t = 0.1;

vec3 nrm(vec2 pos)
{
  float t = min(1.0,length(pos) / R);
  vec3 nP = vec3(0.0,0.0,1.0);          // vertical vector used inside, where the noise is unconstrained
  vec3 nR = vec3(normalize(pos),0.0);   // radial vector constraining the noise along a curve (circle)
  return normalize(nR * t + nP * (1.0-t));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // pixel coordinates
    vec2  p    = (-iResolution.xy + 2.0*fragCoord)/iResolution.y;
    vec3  col  = vec3(0.0);
    float tsrf = length(p);
    vec3  pos  = vec3(p,0.0);
    // ----------
    const float noise_scale = 4.0;
    // amount of displacement
    float d = (0.07 + 0.4 * max(0.0,1.0 - length(p)/R)) * sin(2.0*iTime);
    // number of steps
    const int N = 32;
    // integrate displacement
    float dd = d / float(N);
    vec3  cn;
    vec4  n3  = noised( vec3(pos)*noise_scale );
    vec3  nr  = nrm(pos.xy);
    float tau = n3.x; // iso-contour value
    for (int i=0;i<N;++i) {
      cn  = cross(nr,n3.yzw);
      pos = pos + cn.xyz * dd;
      nr  = nrm(pos.xy);
      n3  = noised( vec3(pos)*noise_scale );
      // project on noise iso-contour
      vec3 g = n3.yzw*noise_scale;
      pos    = pos - g * (n3.x-tau) / (dot(g,g)+1e-30);
    }

    col  = vec3(pattern(pos.xy));

	fragColor = vec4( col, 1.0 );
}