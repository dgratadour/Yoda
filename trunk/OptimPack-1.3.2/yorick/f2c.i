func f2c_warn(s) { write, format="WARNING: %s\n", s; }
func f2c_load(file)
{
  if (structof(file) == string) file = open(file);
  buf = rdline(file, 10);
  while (buf(0)) grow, buf, rdline(file, numberof(buf));
  return buf(where(buf));
}

func f2c_save(filename, buf, overwrite=)
{
  if (! overwrite && open(filename, "r", 1)) error, "file already exists";
  write, open(filename, "w"), format="%s\n", linesize=100000, buf;
}

func f2c(file, trunc=)
/*

   Assume the source lines follows the "Standard Fixed Format":
    - The first 72 columns of each line are scanned.
    - The first five columns must be blank or contain a numeric label. 
    - Continuation lines are identified by a nonblank, nonzero in column 6. 
    - Short lines are padded to 72 characters. 
    - Long lines are truncated.
    - A line with a 'c', 'C', or '*' in column one is a comment line.
    - A totally blank line is a comment line.

 
 */
{
  nil = string(0);
  code = f2c_load(file);

  /* delete exceeding spaces */
  code = regsub("[ \t]+$", code, "");

  /* check for bad characters */
  if (anyof(regmatch("[\t]", code))) error, "some bad characters";
  
  /* Isolate comments from code (a totally blank line or a line with a
     'c', 'C', or '*' in column one is a comment line). */
  is_comment = regmatch("^([cC*]|$)", code);
  comment_lines = where(is_comment);
  comment = code(comment_lines);
  code_lines = where(! is_comment);
  code = code(code_lines);

  /* Convert FORTRAN comments to C-style comments. */
  if (is_array(comment)) {
    comment = strpart(comment, 2:0);
    n = numberof(comment);
    length = strlen(comment);
    (continued = array(int, n))(1:-1) = (comment_lines(dif) == 1);
    new = 1n;
    for (i=1 ; i<=n ; ++i) {
      if (continued(i)) {
        comment(i) = (new ? "/*" : " *") + comment(i);
        new = 0n;
      } else if (length(i)) {
        comment(i) = (new ? "/*" : " *") + comment(i) + " */";
        new = 1n;
      } else {
        comment(i) = (new ? nil : " */");
        new = 1n;
      }
    }
    continued = new = length = [];
  }
  
  /* Deal with lines longer that 72 characters */
  if (max(strlen(code)) > 72) {
    if (trunc) {
      f2c_warn, "truncating lines longer than 72 characters";
      code = strpart(code, 1:72);
    } else {
      f2c_warn, "some lines are longer than 72 characters";
    }
  }

  /* Join continuation lines in code (continuation lines are identified by
     a nonblank, nonzero in column 6). */
  local tail;
  i = where(regmatch("^     [^ 0] *(.*)", code, /* match0 */, tail));
  if (is_array(i)) {
    n = numberof(i);
    if (min(i) <= 1) error, swrite(format="bad continuation line (%d)",
                                   code_lines(min(i)));
    code(i) = nil;
    kp1 = 0;
    for (j=1 ; j<=n ; ++j) {
      if ((k = i(j)) != kp1) l = k - 1;
      code(l) += tail(k);
      kp1 = k + 1;
    }
    i = [];
  }

  /* Locate labels (the first five columns in code lines must be blank or
     contain a numeric label). */
  label_cols = strpart(code, 1:5);
  i = where(regmatch("^ *[0-9]+ *$", label_cols));
  if (is_array(i)) {
    label_value = array(long, numberof(i));
    sread, label_cols(i), label_value;
    label_lines = code_lines(i);
  } else {
    label_value = label_lines = [];
  }
  label_cols = [];
  
  /* Convert code to lower case. */
  code = strlower(code);

  /* Simplify some FORTRAN keywords. */
  code = regsub("([^a-z0-9])go +to([^a-z0-9])", code, "\\1goto\\2", all=1);  
  code = regsub("([^a-z0-9])end *(do|if)$", code, "\\1end");
  

  /* replace some FORTRAN intrinsics */
  code = regsub("^ *subroutine +\(.*\)", code, "void \\1\n{", all=1);
  code = regsub(" *\\.eq\\. *", code, " == ", all=1);
  code = regsub(" *\\.ne\\. *", code, " != ", all=1);
  code = regsub(" *\\.lt\\. *", code, " < ", all=1);
  code = regsub(" *\\.le\\. *", code, " <= ", all=1);
  code = regsub(" *\\.gt\\. *", code, " > ", all=1);
  code = regsub(" *\\.ge\\. *", code, " >= ", all=1);
  code = regsub("\\.true\\.", code, " TRUE ", all=1);
  code = regsub("\\.false\\.", code, " FALSE ", all=1);


  return merge(comment, code, is_comment);
}

/* IF (arithmetic)
     IF ( expr ) s1, s2, s3
     -> { type __tmp = expr;
          if (__tmp < (type)0) goto s1;
          if (__tmp == (type)0) goto s2;
          goto s3;
        }

   IF (block)
     IF ( e1 ) THEN
       ...
     ELSE IF ( e2 ) THEN
       ...
     ELSE
       ...
     END IF
     -> if (e1) {
          ...;
        } else if (e2) {
          ...;
        } else {
          ...;
        }

   
   IF (logical)
     IF ( e ) statement
     -> if (e) statement;


   DO label [,] variable = start, stop [, incr ]

   The expressions e1, e2, and e3 are evaluated. If e3 is not present, its
   value is assumed to be one.

   The iteration count is established as the value of the expression: 
     MAX (INT ((e2 - e1 + e3) / e3 ), 0) 

   The iteration count is zero if either of the following is true: 

     e1 \> e2 and e3 \> zero. 
     e1 < e2 and e3 < zero.

   The iteration count is tested, and, if it is greater than zero, the
   range of the DO loop is executed.

   -> {
         TYPE __stop = e2, __incr = e3;
         for (var = e1 ;
              (__incr > (TYPE)0 ? var <= __stop : var >= __stop) ;
              var += __incr) {
              ...;
         }
      }

*/

#if 0
/*
 * COMMENT
 * CONTINUATION
 * [LABEL] ASSIGN_STMT
 * [LABEL] IF_STMT
 * [LABEL] DO_STMT
 * [LABEL] RETURN_STMT
 */
func f2c_
{
  extern _f2c_re_label;
  local label_values;
  label_cols = strpart(code, 1:5);
  i = where(regmatch("^ *([0-9]+) *$", strpart(code, 1:5), /* match0 */,
                     label_values));
  if (is_array(i)) {
    label_values = array(long, numberof(i));
    sread, label_cols(i), label_values;
    label_lines = code_lines(i);

    label_nrefs = h_new();
    label_lines = h_new();
    n = numberof(i);
    for (j=1 ; j<=n ; ++j) {
      key = swrite(format="%d", label_values(i));
      if (h_has(label_nrefs, key)) error, "duplicate label ("+key+")";
      h_set, label_nrefs, key, 0;
      h_set, label_lines, key, code_lines(i(j));
    }
  } else {
    label_values = label_lines = [];
  }
  
  
  regmatch, ;
}

local _f2c_column, _f2c_line;

/* f2c_parse_line - must be called after dealing with commentary lines and
   continuation lines and after some FORTRAN keywords have been simplified
   ("go to" -> "goto", "end..." -> "end") */
func f2c_parse_line(linestr, linechar, linetype)
{

  column = 0;
  index = bufchar + 1;
  buftype = F2C_TYPE_TABLE(index);
  buf_is_type = F2C_TYPE_TABLE(index);
  
  
  extern _f2c_column, _f2c_linenumber, _f2c_line;
  
  /* skip spaces */
  while ((c = linechar(++_f2c_column)) == ' ') /* noop */ ;

  //if (_f2c_column == 1) {
  //  /* Comment or end of file.  */
  //  if (c == 'c' || c == 'C' || c == '*') {
  //    return COMMENT;
  //  }
  //  if (c == '\n') {
  //    ++_f2c_linenumber;
  //    return NEWLINE;
  //  }
  //  if (c == '\0') {
  //    return EOF;
  //  }
  //  error, swrite(format="illegal character in first column at line %d",
  //                _f2c_linenumber);
  //}
  if (_f2c_column >= 6) {
    i1 = i2 = _f2c_column;
    type = F2C_TYPE(c + 1);
    if (type == F2C_ALPHA) {
      while ((x = linetype(++i2)) >= F2C_ALPHA && x <= F2C_DOLLAR) /* noop */ ;
      range = i1:i2-1;
      token = strpart(line_str, range);
      if (h_has(F2C_KEYWORDS, token)) {
        return F2C_KEYWORD;
      }
      return F2C_VARIABLE;
    }
    if (type == F2C_DIGIT) {
    }
    if (type == F2C_OPEN || type == F2C_CLOSE) {
    }
    if (type == F2C_PERIOD) {
      next_type = ...;
      if (next_type == F2C_ALPHA) {
        while ((x = linetype(++i2)) == F2C_ALPHA) /* noop */ ;
        if (x == F2C_PERIOD) {
          range = i1:i2;
          token = strpart(line_str, range);
          _f2c_column = i2;
          return F2C_INTRINSIC;
        }
      } else if (next_type == F2C_DIGIT) {
        
      }
      range = i1:i2-1;
       
      }
    }
    

    
      do {
        x = linetype(++i2);
      } while (x >= F2C_ALPHA && x <= );
      
      while (linetype(++i2) == char(++_f2c_column)) == ' ') /* noop */ ;
    
    
  } else {
    
  }
  
  if (_f2c_index < _f2c_length) {
    c = linechar(++_f2c_index);
    type = _f2c_type(c + 1);
  }
}

/*
  -1 - illegal character
   0 - end-of-file or end-of-line
   1 - space
   2 - numerical
   3 - alphabetical (can be first letter in a variable name or part of
       variable) 
   4 - operator or punctuation
*/

local F2C_DIGIT;

func f2c_init(extended)
{
  extern _f2c_type;
  extern F2C_SPACE;  F2C_SPACE = 1;
  extern F2C_ALPHA;  F2C_ALPHA = 2;
  extern F2C_DIGIT;  F2C_DIGIT = 3;
  extern F2C_DOLLAR; F2C_DOLLAR = 4;

  extern F2C_OPEN; F2C_OPEN = ;
  extern F2C_; F2C_ = ;
  
  if (F2C_DIGIT != F2C_ALPHA + 1 ||
      F2C_DOLLAR != F2C_ALPHA + 2) error, "assertion failed";
  
  _f2c_type = array(-1, 256);
  _f2c_type('\0' + 1) = 0;
  _f2c_type('\n' + 1) = 0;
  _f2c_type(' ' + 1) = 1;
  _f2c_type(indgen('0'+1:'9'+1)) = F2C_DIGIT;
  _f2c_type(indgen('A'+1:'Z'+1)) = F2C_ALPHA;
  _f2c_type(indgen('a'+1:'z'+1)) = F2C_ALPHA;
  _f2c_type('=' + 1) = F2C_ASSIGN;
  _f2c_type('+' + 1) = F2C_PLUS;
  _f2c_type('-' + 1) = F2C_MINUS;
  _f2c_type('*' + 1) = F2C_ASTERISK;
  _f2c_type('/' + 1) = F2C_SLASH;
  _f2c_type('(' + 1) = F2C_OPEN;
  _f2c_type(')' + 1) = F2C_CLOSE;
  _f2c_type(',' + 1) = F2C_COMMA;
  _f2c_type('.' + 1) = F2C_PERIOD;
  _f2c_type('\'' + 1) = F2C_QUOTE;

  if (extended) {
    _f2c_type('_' + 1) = F2C_ALPHA;
    _f2c_type('$' + 1) = F2C_NAME;
  }
    

  _F2C_IS_ALPHA = array(0n, 256);
  _F2C_IS_ALPHA(indgen('A'+1:'Z'+1)) = 1n;
  _F2C_IS_ALPHA(indgen('a'+1:'z'+1)) = 4;
 
}
#endif
