if(BUILD_DOXYGEN_DOCS)
   configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
   add_custom_command(OUTPUT
      ${CMAKE_CURRENT_BINARY_DIR}/doxygen/html/index.html
      COMMAND ${DOXYGEN_EXECUTABLE}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      VERBATIM)
   # Runs only once
   add_custom_target(apidoc-all ALL
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen/html/index.html)
   # For manual updates
   add_custom_target(apidoc
      COMMAND ${DOXYGEN_EXECUTABLE}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      VERBATIM)
endif()

if(BUILD_MANUAL)
   set(MANUAL_SRCS
      manual/changelogs.tex
      manual/contributors.tex
      manual/cover.tex
      manual/figures/logo_asc.eps
      manual/figures/logo_px200.eps
      manual/figures/note.eps
      manual/figures/tip.eps
      manual/figures/TU_Signet_CMYK.eps
      manual/IEEEtran_v1.13.bst
      manual/installation.tex
      manual/introduction.tex
      manual/keywords.tex
      manual/license.tex
      manual/setup.tex
      manual/viennafvm.bib
      manual/viennafvm.tex)

   set(MANUAL_OBJS
      "${CMAKE_CURRENT_BINARY_DIR}/viennafvm.pdf"
      "${CMAKE_CURRENT_BINARY_DIR}/manual/viennafvm.bbl"
      "${CMAKE_CURRENT_BINARY_DIR}/manual/viennafvm.blg"
      "${CMAKE_CURRENT_BINARY_DIR}/manual/viennafvm.dvi"
      "${CMAKE_CURRENT_BINARY_DIR}/manual/viennafvm.idx"
      "${CMAKE_CURRENT_BINARY_DIR}/manual/viennafvm.log"
      "${CMAKE_CURRENT_BINARY_DIR}/manual/viennafvm.out"
      "${CMAKE_CURRENT_BINARY_DIR}/manual/viennafvm.ps"
      "${CMAKE_CURRENT_BINARY_DIR}/manual/viennafvm.toc")
   foreach(f IN LISTS MANUAL_SRCS)
      configure_file(${f} ${CMAKE_CURRENT_BINARY_DIR}/${f} COPYONLY)
      if(f MATCHES "(.*)\\.tex$")
         list(APPEND MANUAL_OBJS
            "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_MATCH_1}.aux")
      endif()
   endforeach()

   add_custom_command(OUTPUT ${MANUAL_OBJS}
      COMMAND ${LATEX_COMPILER} viennafvm.tex
      COMMAND ${LATEX_COMPILER} viennafvm.tex
      COMMAND ${BIBTEX_COMPILER} viennafvm
      COMMAND ${LATEX_COMPILER} viennafvm.tex
      COMMAND ${BIBTEX_COMPILER} viennafvm
      COMMAND ${LATEX_COMPILER} viennafvm.tex
      COMMAND ${DVIPS_CONVERTER} -Ppdf -G0 -ta4 viennafvm.dvi
      COMMAND ${PS2PDF_CONVERTER} -dPDFSETTINGS=/prepress
         -dCompatibilityLevel=1.3 -dMAxSubsetPct=100 -dSubsetFonts=true
         -dEmbedAllFonts=true -sPAPERSIZE=a4 -dAutoRotatePages=/None
         -dOptimize=true viennafvm.ps ../viennafvm.pdf
      DEPENDS ${MANUAL_SRCS}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/manual
      COMMENT "Generating manual viennafvm.pdf"
      VERBATIM)

   add_custom_target(manual ALL
      DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/viennafvm.pdf")
endif()
