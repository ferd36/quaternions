= LOG
===
Contains a detailed record of the things I experimented with while writing this code.

* 12/24/15:
    - Tested the "fast inverse square root" from Quake for norms, but the precision I could get out of it was too awful.
    I tried this code from the internet:
    ```cpp
    float Q_rsqrt( float number )
    {
    	long i;
    	float x2, y;
    	const float threehalfs = 1.5F;

    	x2 = number * 0.5F;
    	y  = number;
    	i  = * ( long * ) &y;                       // evil floating point bit level hacking
    	i  = 0x5f3759df - ( i >> 1 );               // what the fuck? 0x5fe6eb50c7b537a9 for double
    	y  = * ( float * ) &i;
    	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
    //	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

    	return y;
    }

    ```

* 12/22/15:
    - Just finished installing Travis CI and coveralls. Valgrind is included in the build too.