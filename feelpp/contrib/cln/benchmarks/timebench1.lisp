(progn
  (format t "~%pi")
  (time (setf (long-float-digits) 3322))
  (time pi)
)

(progn
  (format t "~%gamma not yet implemented in CLISP")
)

(progn
  (format t "~%e")
  (time (exp 1L0))
)

(progn
  (format t "~%sqrt(3)")
  (time (sqrt 3L0))
)

(progn
  (format t "~%exp(log(2))")
  (time (exp (log 2L0)))
)

(progn
  (format t "~%log(exp(2))")
  (time (log (exp 2L0)))
)

(progn
  (format t "~%sin(pi/3)")
  (time (sin (/ pi 3)))
)

(progn
  (format t "~%cos(pi/3)")
  (time (cos (/ pi 3)))
)

(progn
  (format t "~%arcsin(sqrt(3)/2)")
  (time (asin (/ (sqrt 3L0) 2)))
)

(progn
  (format t "~%arccos(sqrt(3)/2)")
  (time (acos (/ (sqrt 3L0) 2)))
)

(progn
  (format t "~%sinh(log(2))")
  (time (sinh (log 2L0)))
)

(progn
  (format t "~%cosh(log(2))")
  (time (cosh (log 2L0)))
)

(progn
  (format t "~%arsinh(pi)")
  (time (asinh pi))
)

(progn
  (format t "~%arcosh(pi)")
  (time (acosh pi))
)
