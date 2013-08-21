;; -*- emacs-lisp -*-
(
  (nil .
    (
      (tab-width . 4)
      (indent-tabs-mode . nil)
      (fill-column . 80)
      (user-company . "Feel++ Consortium")
      (user-mail-address . "christophe.prudhomme@feelpp.org")
      (debian-changelog-mailing-address . "christophe.prudhomme@feelpp.org")
      ))
  (c++-mode .
    (
      (subdirs . nil)
      (indent-tabs-mode . nil)
      (tab-width . 4)
      (show-trailing-whitespace . t)
      (indicate-empty-lines . t)
      (c-basic-offset . 4)
      (c-set-offset 'innamespace  0)
      (add-hook 'write-file-functions 'delete-trailing-whitespace)
      ))
  (c-mode .
    (
      (subdirs . nil)
      (indent-tabs-mode . nil)
      (tab-width . 4)
      (show-trailing-whitespace . t)
      (indicate-empty-lines . t)
      (c-basic-offset . 4)
      (c-set-offset 'innamespace  0)
      (add-hook 'write-file-functions 'delete-trailing-whitespace)
      ))
  )
