
# STYLE GUIDE

Ref: [https://www.fortran90.org/src/best-practices.html](https://www.fortran90.org/src/best-practices.html)
                            
Here is a style guide that we like and that seems to be prevalent 
in most scientific codes (as well as the Fortran standard library).

### Naming Convention:

- Use lowercase for all Fortran constructs (e.g. `do`, `subroutine`, `module`, `type`, `omp parallel` ...).

- Follow short mathematical notation for mathematical variables/functions (e.g. delta_t, Ylm, Gamma, gamma, Enl, Rnl, ...).

- For other names use all lowercase.

- Try to keep names to one or two syllables; if more are required, use underscores to clarify. No camelCase.
Example: sortpair, whitechar, meshexp, numstrings, linspace, meshgrid, argsort, spline, spline_interp, spline_interpolate,
stoperr, stop_error, meshexp_der. 

> Note, for example, “spline interpolation” can be shortened to spline_interpolation, spline_interpolate, spline_interp, spline, but
not to splineint (“int” could mean integration, integer, etc. — too much ambiguity, even in the clear context of a computational code).
This is in contrast to get_argument() where getarg() is perfectly clean and clear.

The above are general guidelines. In general, choosing the right name certainly 
depends on the word being truncated as to whether the first syllable is sufficient. 
Usually it is but clearly not always. Thus some thought should go into step “try to
keep names to 2 syllables or less” since it can really affect the indicativeness and
simplicity. Simple consistent naming rules are a real help in this regard – for both
collaboration and for one’s own sanity when going back to some old code you haven’t 
seen in while.

### Indentation:

- Use 4 spaces indentation (set your source-code editor to put 4 Spaces when you press Tab).
- Code blocks inside every construct (e.g. `do`...`end do`, `module`...`end module`, `subroutine`....`end subroutine`) should be indented.
- One blank line should be in between two code blocks
-  Leave one space between consecutive terms/keywords. E.g. `i  = 3*x + y`, `real :: y`
```fortran
    module types
        integer :: i
        
        do i = 1,10
            i = i
            .......
        end do
    end module    
```

### Columns
If a line stretches more than 120 columns, break it.

### Comments
Inline comments should not start with capital letter. Full comment lines should start with caps. Omp sentinels however should be lowercase: `!$omp`.

### Array declaration
Use `dimension()` format. 
