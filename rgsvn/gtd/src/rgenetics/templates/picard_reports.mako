<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
  <head>
    <meta http-equiv="Content-Type" content="text/xhtml; charset=utf-8" />
    <meta name="generator" content="Galaxy ${program_name} tool output - see http://g2.trac.bx.psu.edu/" />
    <link rel="stylesheet" href="/static/style/base.css" type="text/css" />
  </head>
  <body>
    <div class="document">
      <h3><a href="http://rgenetics.org">Rgenetics</a> tool ${program_name} run at ${timestamp}</h3>
      <b>Your job produced the following outputs - check here for a record of what was done and any unexpected events</b>
      <hr />
      % if file_info:
      <div>
	<b>Output files.</b>
	<table>
	  % for file_url, file_name in file_info:
	  <tr><td><a href="${file_url}">${file_name}</a></td></tr>
	  % endfor
	</table>
      </div>
      <hr />
      % endif
      <div>
	<b>Log of activity</b>
	<hr/>
	${log_data}
      </div>
      <div>
	<p>Note: The freely available 
	  <a href="http://picard.sourceforge.net/command-line-overview.shtml">Picard software</a>
	  generated all outputs reported here. These third party tools were orchestrated by the Galaxy 
	  ${program_name} wrapper and this command line from the Galaxy form:
	</p>
	<hr />
	<div>${command_string}</div>
      </div>
  </body>
</html>
