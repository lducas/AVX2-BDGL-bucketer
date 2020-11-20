#include <iostream>
#include <bitset>
#include "math.h"
#include "fht_lsh.h"
#include <chrono>

int N = 49;

void sample_spherical(int N, float *tab)
{
  float norm = 0;
  float corner_case = (rand()) / static_cast <float> (RAND_MAX);
  if (corner_case < 0.001)
  {
    for (int i = 0; i < N; ++i) tab[i] = 0;
    size_t j =  (static_cast<size_t> (rand())) % N;
    tab[j] = 1;
    return;
  }

  for (int i = 0; i < N; ++i) 
  {
    tab[i] =  (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - .5);
    tab[i] += (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - .5);
    tab[i] += (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - .5);
    tab[i] += (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - .5);
    tab[i] += (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - .5);
    tab[i] += (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - .5);
    norm += tab[i]*tab[i];
  }
  norm = 1/sqrt(norm);

//  tab[0] *= 10;
  for (int i = 0; i < N; ++i) 
  {
    tab[i] *= -norm;
  }
}

float ip(int N, float *tab, float *tab0)
{
  double res = 0;
  for (int i = 0; i < N; ++i) 
  {
    res += tab[i]*tab0[i];
  }
  return res;
}

int main() {

// {
//     float v[70] = {-1.03428, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//     ProductLSH lsh(70, 3, 8112, 4, 0, 1);
//     int32_t h[4];
//     lsh.hash(v, h);
//     if( std::abs(h[0])-1 >= 8112 )
//         std::cerr << h[0] << std::endl;
//     assert( std::abs(h[0])-1 < 8112 );
// }
// return 0;

  // {
  //   // benchmark
  //   int64_t vecs = 10000/10;
  //   int64_t samples = 15887169/100; // 0.95 * 3.2 * (4/3.)**60 / 6
  //   int N = 120;
  //   float v[vecs][N];
  //   for( int i = 0; i < vecs; i++ )
  //       sample_spherical(N, &v[i][0]);
  //   std::cerr << "Construct code" << std::endl;
  //   for( int block = 1; block <= 3; block++ ) {
  //       for( int M = 1; M < 20; M += 5 ) {    
  //           ProductLSH lsh(N, block, 1448, M, 1);
  //           int32_t h[M];
  //           auto time_start = std::chrono::steady_clock::now();
  //           int sm = 0;
  //           for( int i = 0; i < samples; i++ ) {
  //               lsh.hash(&v[i%vecs][0], h);
  //               sm += h[0];
  //           }
  //           std::cerr << block << " " << M << " " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - time_start).count() << " " << sm << std::endl;
    
  //       }
  //   }
  // }
  // return 0;
  for (int N = 64; N < 128; N+=1)
  {

    for (int M = 1; M < 20; M+=15)
    {
      printf("\n\n");
      for (int X = 1; X < 4; X++)
      {

        int64_t samples = (1 << 15);
        int64_t buckets = M * pow(2,  .2015 *N);
        int64_t references = (1 << 14);

        ProductLSH lsh(N, X, buckets, M, 1);
        // actually consider the true number of buckets
        buckets = lsh.codesize;

        std::vector<std::vector<float> > ref(references, std::vector<float>(N));
        int32_t h[M]; 
        std::vector<std::vector<size_t> > buck(buckets);

        float v[N] = {0};

        double ip_pass = 0, ip_all = 0;
        int pass = 0;

        for (int i = 0; i < references; ++i)
        {
          sample_spherical(N, ref[i].data());


          lsh.hash(ref[i].data(), h);

          for (int j = 0; j < M; ++j)
          {
            int hh = h[j];
            assert(hh >= 0);
            assert(hh < buckets);

            buck[hh].push_back(i);
  
          }
        }

        for (int r = 0; r < samples; ++r)
        {
          sample_spherical(N, v);
          lsh.hash(v, h);

          float ipp = ip(N, v, ref[0].data());
          ip_all += ipp*ipp;

          for (int j = 0; j < M; ++j)
          {
            int hh = h[j];
            assert(hh >= 0);
            assert(hh < buckets);
            for (int k : buck[hh])
            {
              float ipp = ip(N, v, ref[k].data());
              ip_pass += ipp*ipp;
              pass ++;              
            }
          }
        }

        printf("\n dim %2d buckets %6d blocks %d multi %d \n", N,(int) buckets, X, M);
        printf("collisions ratio : %6d / %e = %f *expected \n", pass, (1.*references)*samples, (1.*buckets*pass) / (references*samples * M * M));
        printf("std-ip : collision %.4f  all %.4f \n", sqrt(ip_pass/pass) ,sqrt(ip_all/(samples)));
      }
    }
  }
  return 0;
}
