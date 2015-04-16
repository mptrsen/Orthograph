push:
	sed -i -e "s/'[^']\+'/'$(git describe --abbrev=0)'/" Orthograph/Version.pm
	git commit -a
	git log --pretty --date=short > ChangeLog
	git commit --all --amend --no-edit
	git push
