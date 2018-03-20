all:
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/joss-logo.png
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/latex.template
	pandoc \
    -V repository="https://github.com/pynucastro/pynucastro" \
    -V archive_doi="http://dx.doi.org/10.5281/zenodo.1202434" \
    -V paper_url="http://joss.theoj.org/papers/10.21105/joss.00588" \
    -V formatted_doi="10.21105/joss.00588" \
    -V review_issue_url="https://github.com/openjournals/joss-reviews/issues/588" \
    -V graphics="true" \
	-V journal_name="Journal of Open Source Software" \
    -V issue="23" \
    -V volume="3" \
    -V page="588" \
    -V logo_path="joss-logo.png" \
    -V year="2018" \
    -V submitted="15 February 2018" \
    -V published="18 March 2018" \
    -V citation_author="Willcox et al." \
    -V paper_title="pynucastro: an interface to nuclear reaction rates and code generator for reaction network equations" \
    paper.md --template latex.template --filter pandoc-citeproc --bibliography paper.bib -o paper.pdf -V geometry:margin=1in

tex:
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/joss-logo.png
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/latex.template
	pandoc \
    -V repository="https://github.com/pynucastro/pynucastro" \
    -V archive_doi="http://dx.doi.org/10.5281/zenodo.1202434" \
    -V paper_url="http://joss.theoj.org/papers/10.21105/joss.00588" \
    -V formatted_doi="10.21105/joss.00588" \
    -V review_issue_url="https://github.com/openjournals/joss-reviews/issues/588" \
    -V graphics="true" \
	-V journal_name="Journal of Open Source Software" \
    -V issue="23" \
    -V volume="3" \
    -V page="588" \
    -V logo_path="joss-logo.png" \
    -V year="2018" \
    -V submitted="15 February 2018" \
    -V published="18 March 2018" \
    -V citation_author="Willcox et al." \
    -V paper_title="pynucastro: an interface to nuclear reaction rates and code generator for reaction network equations" \
    paper.md --template latex.template --filter pandoc-citeproc --bibliography paper.bib -o paper.tex -V geometry:margin=1in
