/*
  Accompanying code, Improving Curl Noise, SIGGRAPH Asia 2025
  ---

  This demonstrates curl noise along a torus surface, obtained as
  the cross product of a solid noise gradient and the gradient of
  the torus distance field.

  Accuracy is greatly improved by using re-projection onto the
  iso-contours of the noise and surface at every step.

*/

// enables/disables projection
#define USE_PROJECTION 1

// noise scale (higher => more swirls)
const float noise_scale = 6.0;

// number of steps
const int N = 128; // 32, 128, 2048 (reference)

// -------------------------------------------------------
// -------------------------------------------------------
// Support code for noise and torus rendering, skip over
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

// routine for ray-torus intersection, https://www.shadertoy.com/view/4sBGDy
// The MIT License
// Copyright © 2014 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
float iTorus( in vec3 ro, in vec3 rd, in vec2 torus )
{
    float Ra2 = torus.x*torus.x;
    float ra2 = torus.y*torus.y;
    float m = dot(ro,ro);
    float n = dot(ro,rd);
    float k = (m - ra2 - Ra2)/2.0;
    float a = n;
    float b = n*n + Ra2*rd.z*rd.z + k;
    float c = k*n + Ra2*ro.z*rd.z;
    float d = k*k + Ra2*ro.z*ro.z - Ra2*ra2;
    //----------------------------------
    float p = -3.0*a*a     + 2.0*b;
    float q =  2.0*a*a*a   - 2.0*a*b   + 2.0*c;
    float r = -3.0*a*a*a*a + 4.0*a*a*b - 8.0*a*c + 4.0*d;
    p /= 3.0;
    r /= 3.0;
    float Q = p*p + r;
    float R = 3.0*r*p - p*p*p - q*q;

    float h = R*R - Q*Q*Q;
    float z = 0.0;
    if( h < 0.0 )
    {
        float sQ = sqrt(Q);
        z = 2.0*sQ*cos( acos(R/(sQ*Q)) / 3.0 );
    }
    else
    {
        float sQ = pow( sqrt(h) + abs(R), 1.0/3.0 );
        z = sign(R)*abs( sQ + Q/sQ );
    }
    z = p - z;
    //----------------------------------
    float d1 = z   - 3.0*p;
    float d2 = z*z - 3.0*r;
    if( abs(d1)<1.0e-5 )
    {
        if( d2<0.0 ) return -1.0;
        d2 = sqrt(d2);
    }
    else
    {
        if( d1<0.0 ) return -1.0;
        d1 = sqrt( d1/2.0 );
        d2 = q / d1;
    }
    //----------------------------------
    float result = 1e20;
    h = d1*d1 - z + d2;
    if( h>0.0 )
    {
        h = sqrt(h);
        float t1 = -d1 - h - a;
        float t2 = -d1 + h - a;
             if( t1>0.0 ) result=t1;
        else if( t2>0.0 ) result=t2;
    }
    h = d1*d1 - z - d2;
    if( h>0.0 )
    {
        h = sqrt(h);
        float t1 = d1 - h - a;
        float t2 = d1 + h - a;
             if( t1>0.0 ) result=min(result,t1);
        else if( t2>0.0 ) result=min(result,t2);
    }

    return result;
}
// df(x)/dx
vec3 nTorus( in vec3 pos, vec2 tor )
{
    return normalize( pos*(dot(pos,pos)- tor.y*tor.y - tor.x*tor.x*vec3(1.0,1.0,-1.0)));
}
// distance field
float dTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xy)-t.x,p.z);
  return length(q)-t.y;
}
// -------------------------------------------------------
// -------------------------------------------------------
// end of support code
// -------------------------------------------------------
// -------------------------------------------------------

// generates a simple square pattern
float pattern(vec3 p)
{
  vec3 f = fract(p*5.0);
  vec3 d = abs(0.5-f);
  return 3.0*min(d.x,min(d.y,d.z));
}

// shader entry point
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // camera movement
    //float an = 1.4;
    float an = 1.4 + sin(iTime * 0.4);
    vec3 ro  = vec3( 0.0, -1.9 + 1.5*cos(an), 1.3);
    vec3 ta  = vec3( 0.0, 0.0, -0.6);
    // camera matrix
    vec3 ww  = normalize( ta - ro );
    vec3 uu  = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
    vec3 vv  = normalize( cross(uu,ww));
    vec3 tot = vec3(0.0);
    // pixel coordinates
    vec2 p   = (-iResolution.xy + 2.0*fragCoord)/iResolution.y;
    // create view ray
    vec3 rd  = normalize( p.x*uu + p.y*vv + 1.5*ww );
    // raytrace torus
    vec2 torus = vec2(1.0,0.45 /* + 0.04 * (1.0 - cos(iTime)) */ );
    float t    = iTorus( ro, rd, torus );
    // shading/lighting
    vec3 col = vec3(0.0);
    if( t>0.0 && t<100.0 ) {
        // there is a hit, get the normal
        vec3 pos = ro + t*rd;
        vec3 nor = nTorus( pos, torus );
        // shade
        float dif = clamp( dot(nor,vec3(0.57703)), 0.0, 1.0 );
        float amb = 0.5;
        // ----------
        // advection along curl-noise
        // amount of displacement
        float d = cos(iTime * 0.5) * 0.35;
        float dd = d / float(N);
        vec3  cn;
        vec4  n4  = noised( pos*noise_scale );
        nor       = nTorus( pos, torus );
        float tau = n4.x; // iso-contour value
        for (int i=0;i<N;++i) {
          cn  = cross(nor,n4.yzw);
          pos = pos + cn * dd;
          n4  = noised( pos*noise_scale );
          nor = nTorus( pos, torus );
#if USE_PROJECTION==1
          vec3 g = n4.yzw*noise_scale;
          pos = pos - g * (n4.x-tau) / (dot(g,g)+1e-30);  // project back along the noise iso-contour
          pos = pos - nor * (dTorus(pos,torus)); // project back along the torus surface (iso-contour)
#endif
        }
        vec3 clr = vec3(pattern(pos));
        // ----------
        // done
        col       = vec3((dif+amb) * clr);
    } else {
        // no hit, background color
        col = vec3(0.9);
    }
    fragColor = vec4( sqrt(col), 1.0 );
}
