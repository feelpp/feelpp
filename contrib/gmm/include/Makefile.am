
includedir=${prefix}/include

EXTRA_HEADERS=\
          dal_algobase.h\
          dal_basic.h\
          dal_except.h\
          dal_ref.h\
          dal_std.h\
          gmm_blas.h\
          gmm_blas_interface.h\
          gmm_condition_number.h\
          gmm_conjugated.h\
          gmm_def.h\
          gmm_dense_Householder.h\
          gmm_dense_lu.h\
          gmm_dense_qr.h\
          gmm_dense_sylvester.h\
          gmm_domain_decomp.h\
          gmm.h\
          gmm_inoutput.h\
          gmm_interface_bgeot.h\
          gmm_interface.h\
          gmm_iter.h\
          gmm_iter_solvers.h\
          gmm_kernel.h\
          gmm_lapack_interface.h\
          gmm_matrix.h\
          gmm_modified_gram_schmidt.h\
          gmm_opt.h\
          gmm_precond_diagonal.h\
          gmm_precond.h\
          gmm_precond_ildlt.h\
          gmm_precond_ildltt.h\
          gmm_precond_ilu.h\
          gmm_precond_ilut.h\
          gmm_precond_ilutp.h\
          gmm_precond_mr_approx_inverse.h\
          gmm_real_part.h\
          gmm_scaled.h\
          gmm_solver_bfgs.h\
          gmm_solver_bicgstab.h\
          gmm_solver_cg.h\
          gmm_solver_constrained_cg.h\
          gmm_solver_gmres.h\
          gmm_solver_idgmres.h\
          gmm_solver_qmr.h\
          gmm_solver_Schwarz_additive.h\
          gmm_sub_index.h\
          gmm_sub_matrix.h\
          gmm_sub_vector.h\
          gmm_superlu_interface.h\
          gmm_transposed.h\
          gmm_tri_solve.h\
          gmm_vector.h\
          gmm_vector_to_matrix.h

EXTRA_DIST=$(EXTRA_HEADERS)

install-data-hook:
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(includedir)
	@list='$(EXTRA_HEADERS)'; for p in $$list; do \
	  if test -f "$$p"; then d=; else d="$(srcdir)/"; fi; \
	  f="`echo $$p | sed -e 's|^.*/||'`"; \
	  echo " $(INSTALL_DATA) $$d$$p $(DESTDIR)$(includedir)/$$f"; \
	  $(INSTALL_DATA) $$d$$p $(DESTDIR)$(includedir)/$$f; \
	done

uninstall-hook:
	@$(NORMAL_UNINSTALL)
	@list='$(EXTRA_HEADERS)'; for p in $$list; do \
	  f="`echo $$p | sed -e 's|^.*/||'`"; \
	  echo " rm -f $(DESTDIR)$(includedir)/$$f"; \
	  rm -f $(DESTDIR)$(includedir)/$$f; \
	done

