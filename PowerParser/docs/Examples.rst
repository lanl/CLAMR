========
Examples
========

---------
Example 1
---------

This example will go through writing a quadratic equation subroutine in the input file. One
reason for doing this is to use more complicated examples to help find parser bugs before
the users find them. The second reason is to provide a step by step example of writing a
subroutine including the problems encountered.

The quadratic equation is:

.. math::

   ax^2 + bx + c = 0

and its solution is:

.. math::

   x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}

We start the quadratic subroutine with the following line::

   subroutine quad_eq($a, $b, $c, $root1, $root2)

Next, it is good practice to initialize the output variables.
::

      $root1 = 0.
      $root2 = 0.

Then we define the argument to the sqrt function and make sure it is greater than or equal
to 0. It is good practice to make checks like this in the user input file so as to get more
informative error messages and make it easier to fix the input file.
::

      $rarg = ($b**2 - 4*$a*$c)
      if ($rarg .lt. 0.) then
         fatal_error sub quad_eq, sqrt argument < 0.
      endif

Now define the numerator and denominator for the roots and make sure that the denominator
is not 0.
::

      $num = (-$b + sqrt($rarg))
      $dem = (2 * $a)
      if ($dem .eq. 0.) then
         fatal_error sub quad_eq, denominator (2a) = 0.
      endif

Finally calculate the roots.
::

      $root1 = $num / $dem
      $num = (-$b - sqrt($rarg))
      $root2 = $num / $dem
   end subroutine

Notice that a mistake has been made in calculating the roots. The parentheses around the
math expressions are missing. This is one of the more common errors in developing input
files. Another common error is missing dollar signs in front of variable names. The corrected
input is::

      $root1 = ($num / $dem)
      $num = (-$b - sqrt($rarg))
      $root2 = ($num / $dem)
   end subroutine

The final subroutine is given here::

   !===============================================================================
   ! Subroutine to solve the quadratic equation.
   !===============================================================================
   subroutine quad_eq($a, $b, $c, $root1, $root2)
      $root1 = 0.
      $root2 = 0.
      $rarg = ($b**2 - 4*$a*$c)
      if ($rarg .lt. 0.) then
         fatal_error sub quad_eq, sqrt argument < 0.
      endif
      $num = (-$b + sqrt($rarg))
      $dem = (2 * $a)
      if ($dem .eq. 0.) then
         fatal_error sub quad_eq, denominator (2a) = 0.
      endif
      $root1 = ($num / $dem)
      $num = (-$b - sqrt($rarg))
      $root2 = ($num / $dem)
   end subroutine

The first attempt at calling the quadratic subroutine was::

   call quad_eq(2, 4, -30, $x1, $x2)

This had 2 problems. The first is that $x1 and $x2 have not been defined. It is natural to
want to use them in this way since we know that the subroutine is going to define them. But
the current parser will not allow that. This perhaps will be fixed sometime. The corrected
input is:

   $x1 = 0
   $x2 = 0
   call quad_eq(2, 4, -30, $x1, $x2)

The second problem is that there was a bug in the parser that the unary minus sign in front
of the number 30 is not being applied correctly. This has been fixed. Unfortunately problems
like this might crop up as the parser is being debugged. Users are encouraged to report the
bug and to attempt a workaround. In this case, treating the -30 as a math expression and
enclosing it in parentheses worked, i.e.
::

   $x1 = 0
   $x2 = 0
   call quad_eq(2, 4, (-30), $x1, $x2)

It would have also worked to define a new variable equal to -30 and use the variable in the
subroutine call.

For this example, the roots, $x1 and $x2, end up being 3 and -5.

It is good practice to verify that the fatal_error commands work as expected. The division
by zero can be checked by calling the subroutine with $a = 0
::

   call quad_eq(0, 4, -30, $x1, $x2)

The sqrt of a negative number error can be checked by changing -30 to +30.
::

   call quad_eq(2, 4, 30, $x1, $x2)

