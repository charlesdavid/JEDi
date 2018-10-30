package jedi;

import java.io.IOException;
import java.io.StringReader;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.util.Stack;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * The class Fortran Format
 * 
 * FortranFormat Version 1.0, written by Kevin J. Theisen
 * 
 * Copyright (c) 2009 iChemLabs, LLC. All rights reserved.
 * 
 * $Revision: 793 $ $Author: kevin $ $LastChangedDate: 2009-03-14 20:03:16 -0400 (Sat, 14 Mar 2009) $
 * 
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 2. Redistributions in binary form must reproduce
 * the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 3. Neither the name
 * of the iChemLabs nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

public class FortranFormat
	{
		public static Vector<Object> obj;

		/**
		 * Read function similar to Fortran implementation.
		 * 
		 * @param data   is the data to be parsed
		 * @param format is the format specification
		 * 
		 * @return a Vector<object> of all the parsed data as Java objects
		 * 
		 * @throws ParseException the parse exception
		 * @throws IOException    Signals that an I/O exception has occurred.
		 */

		public static Vector<Object> read(String data, String format)
			{
				FortranFormat ff = new FortranFormat(format);
				obj = ff.parse(data);
				return obj;
			}

		/**
		 * Write function similar to the Fortran implementation.
		 * 
		 * @param objects is the vector of objects to be formatted
		 * @param format  is the format specification
		 * 
		 * @return the formatted string
		 * 
		 * @throws ParseException the parse exception
		 * @throws IOException    Signals that an I/O exception has occurred.
		 */

		public static String write(Vector<Object> objects, String format)
			{
				FortranFormat ff = new FortranFormat(format);
				return ff.format(objects);
			}

		/**
		 * The Class Unit. Holds a single Edit Descriptor.
		 */
		private class Unit
			{

				/** The type. */
				private String type;

				/** The length 'w'. */
				private int length;
				/** The decimal length 'd'. */
				private int decimalLength;
				/** The exponent length 'e'. */
				private int exponentLength;

				/**
				 * Instantiates a new unit.
				 * 
				 * @param type   the type
				 * @param length the length 'w'
				 */
				public Unit(String type, int length)
					{
						this.type = type;
						this.length = length;
					}

				/*
				 * (non-Javadoc)
				 * 
				 * @see java.lang.Object#toString()
				 */
				@Override
				public String toString()
					{
						return type + length + (decimalLength > 0 ? "." + decimalLength : "") + (exponentLength > 0 ? "E" + exponentLength : "") + " ";
					}

			}

		/** A constant for the Repeatable Integer Edit Descriptor. */
		protected static String R_INTEGER = "I";

		/** A constant for the Repeatable Real Decimal Edit Descriptor. */
		protected static String R_REAL_DECIMAL = "F";

		/** A constant for the Repeatable Real Exponent Edit Descriptor. */
		protected static String R_REAL_EXPONENT = "E";

		/** A constant for the Repeatable Real Scientific Edit Descriptor. */
		protected static String R_REAL_SCIENTIFIC = "ES";

		/** A constant for the Repeatable Real Engineering Edit Descriptor. */
		protected static String R_REAL_ENGINEERING = "EN";

		/** A constant for the Repeatable Real Double Edit Descriptor. */
		protected static String R_REAL_DOUBLE = "D";

		/** A constant for the Repeatable Real Decimal (Redundant) Edit Descriptor. */
		protected static String R_REAL_DECIMAL_REDUNDANT = "G";

		/** A constant for the Repeatable Logical Edit Descriptor. */
		protected static String R_LOGICAL = "L";

		/** A constant for the Repeatable Character Edit Descriptor. */
		protected static String R_CHARACTER = "A";

		/** A constant for the Non-repeatable Horizontal Positioning Edit Descriptor. */
		protected static String NR_POSITIONING_HORIZONTAL = "X";

		/** A constant for the Non-repeatable Positioning Tab Edit Descriptor. */
		protected static String NR_POSITIONING_TAB = "T";

		/** A constant for the Non-repeatable Positioning Tab Left Edit Descriptor. */
		protected static String NR_POSITIONING_TAB_LEFT = "TL";

		/** A constant for the Non-repeatable Positioning Tab Right Edit Descriptor. */
		protected static String NR_POSITIONING_TAB_RIGHT = "TR";

		/** A constant for the Non-repeatable Vertical Positioning Edit Descriptor. */
		protected static String NR_POSITIONING_VERTICAL = "/";

		/** A constant for the Non-repeatable Scanning Control Edit Descriptor. */
		protected static String NR_FORMAT_SCANNING_CONTROL = ":";

		/** A constant for the Non-repeatable Sign Control Default Tab Edit Descriptor. */
		protected static String NR_SIGN_CONTROL_COMPILER = "S";

		/** A constant for the Non-repeatable Sign Control Always Use Pluses Tab Edit Descriptor. */
		protected static String NR_SIGN_CONTROL_POSITIVE_ALWAYS = "SP";

		/** A constant for the Non-repeatable Sign Control Never Use Pluses Tab Edit Descriptor. */
		protected static String NR_SIGN_CONTROL_POSITIVE_NEVER = "SS";

		/** A constant for the Non-repeatable Remove Blanks Edit Descriptor. */
		protected static String NR_BLANK_CONTROL_REMOVE = "BN";

		/** A constant for the Non-repeatable Blanks-as-Zeros Edit Descriptor. */
		protected static String NR_BLANK_CONTROL_ZEROS = "BZ";

		/** The format specification. */
		private String format;

		/** The parsed Edit Descriptors. */
		private Vector<Unit> units;

		/** The char to use when skipping positions during write. */
		private char positioningChar = ' ';

		/** Use this to set whether or not to append a return line to the end of the generated string during write. */
		private boolean addReturn = false;

		/**
		 * Instantiates a new FortranFormat object.
		 * 
		 * @param format is the format specification string
		 * 
		 * @throws ParseException the parse exception
		 */
		public FortranFormat(String format)
			{
				this.format = format;
				parseFormat();
			}

		/**
		 * Parses the format specification string.
		 * 
		 * @throws ParseException the parse exception
		 */
		private void parseFormat()
			{
				String noParen = removeParenthesis(format);
				final StringTokenizer st = new StringTokenizer(noParen, " ,");
				units = new Vector<Unit>(st.countTokens());
				while (st.hasMoreTokens())
					{
						final String s = st.nextToken();
						boolean reachedType = false, hasDecimal = false, hasExponent = false;
						final StringBuffer before = new StringBuffer(), type = new StringBuffer(), decimal = new StringBuffer(), exponent = new StringBuffer();
						StringBuffer after = new StringBuffer();
						for (int i = 0; i < s.length(); i++)
							{
								if (s.charAt(i) == '.')
									{
										hasDecimal = true;
									} else if (reachedType && s.charAt(i) == 'E')
									{
										hasExponent = true;
									} else if (Character.isLetter(s.charAt(i)) || s.charAt(i) == '/')
									{
										type.append(s.charAt(i));
										reachedType = true;
									} else
									{
										if (hasExponent)
											{
												exponent.append(s.charAt(i));
											} else if (hasDecimal)
											{
												decimal.append(s.charAt(i));
											} else if (reachedType)
											{
												after.append(s.charAt(i));
											} else
											{
												before.append(s.charAt(i));
											}
									}
							}
						int repeats = before.length() == 0 ? 1 : Integer.parseInt(before.toString());
						if (type.toString().equals(NR_POSITIONING_HORIZONTAL))
							{
								after = before;
								repeats = 1;
							}
						if (type.toString().equals("E"))
							{
								if (exponent.length() == 0)
									{
										exponent.append('2');
									}
							}
						for (int i = 0; i < repeats; i++)
							{
								Unit u = new Unit(type.toString(), after.length() == 0 ? 0 : Integer.parseInt(after.toString()));
								if (decimal.length() != 0)
									{
										u.decimalLength = Integer.parseInt(decimal.toString());
									}
								if (exponent.length() != 0)
									{
										u.exponentLength = Integer.parseInt(exponent.toString());
									}
								units.add(u);
							}
					}
			}

		/**
		 * Adds the commas to the correct places.
		 * 
		 * @param withoutCommas is the string to be analyzed
		 * 
		 * @return a string with proper comma placement
		 */
		private static String addCommas(String withoutCommas)
			{
				StringBuffer sb = new StringBuffer();
				boolean hitE = false;
				boolean lastWasChar = true;
				boolean foundNotNum = false;
				for (int i = 0; i < withoutCommas.length(); i++)
					{
						char c = withoutCommas.charAt(i);
						if (c == '.' || Character.isDigit(c))
							{
								sb.append(c);
								lastWasChar = false;
								if (i != 0 && withoutCommas.charAt(i - 1) == ',')
									{
										foundNotNum = false;
									}
							} else
							{
								if (foundNotNum && !lastWasChar && c != ',' && c != '(' && c != ')' && c != 'X' && i != 0 && withoutCommas.charAt(i - 1) != ','
										&& !(c == 'E' && hitE))
									{
										sb.append(',');
										hitE = false;
									}
								if (c == 'E')
									{
										hitE = true;
									}
								foundNotNum = true;
								lastWasChar = true;
								sb.append(c + (c == '/' ? "," : ""));
							}
					}
				return sb.toString();
			}

		/**
		 * Removes the parenthesis from the specification string.
		 * 
		 * @param withParen is the string with parenthesis
		 * 
		 * @return an equivalent string without any parenthesis
		 * 
		 * @throws ParseException the parse exception
		 */
		private static String removeParenthesis(String withParen)
			{
				StringBuffer sb = new StringBuffer();
				int open = withParen.indexOf('(');
				try
					{
						if (open == -1)
							{
								throw new ParseException(
										"Fortran format specification strings must begin with an open parenthesis '(' and end with a close parenthesis ')'. Blank spaces are tolerated before an open parenthesis and any characters are tolerated after a close parenthesis. No characters outside of the root parenthesis affect the format specification.",
										0);
							}
						int close = findClosingParenthesis(withParen, open);
						String before = withParen.substring(0, open);
						String content = withParen.substring(open + 1, close);
						for (int i = 0; i < before.length(); i++)
							{
								if (before.charAt(i) != ' ')
									{
										throw new ParseException("Only spaces may precede the root parenthesis.", 0);
									}
							}
						content = content.replaceAll(" ", "");
						sb.append(recursiveRemoveParenthesis(content));

					} catch (ParseException e)
					{
						System.out.println("ParseException thrown. PDB file not compatible with chosen format");
						e.getMessage();
						e.getStackTrace();
					}

				return sb.toString();
			}

		/**
		 * Recursively remove parenthesis. Also forwards output to addCommas().
		 * 
		 * @param content is the string that may contain parenthesis
		 * 
		 * @return an equivalent string without any parenthesis
		 * 
		 * @throws ParseException the parse exception
		 */
		private static String recursiveRemoveParenthesis(String content)
			{
				StringBuffer sb = new StringBuffer();
				StringBuffer num = new StringBuffer();
				for (int i = 0; i < content.length(); i++)
					{
						char c = content.charAt(i);
						if (c == '(')
							{
								int closing = findClosingParenthesis(content, i);
								int repeat = 1;
								if (num != null && num.length() > 0)
									{
										String s = num.toString();
										if (s.startsWith(","))
											{
												s = s.substring(1);
											}
										boolean isNum = true;
										for (int j = 0; j < s.length(); j++)
											{
												if (!Character.isDigit(s.charAt(j)))
													{
														isNum = false;
													}
											}
										if (isNum)
											{
												repeat = Integer.parseInt(s);
											} else
											{
												sb.append(s);
											}
									}
								String total = recursiveRemoveParenthesis(content.substring(i + 1, closing));
								i += closing - i;
								for (int j = 0; j < repeat; j++)
									{
										if (sb.length() > 0)
											{
												sb.append(',');
											}
										sb.append(total);
									}
								num = new StringBuffer();
							} else
							{
								num.append(c);
							}
					}
				if (num.length() > 0 && sb.length() > 0)
					{
						sb.append(',');
					}
				sb.append(num);
				String complete = sb.toString();
				if (complete.indexOf(',') == -1 || complete.indexOf('/') != -1)
					{
						complete = addCommas(complete);
					}
				return complete;
			}

		/**
		 * Find the closing parenthesis to a given open parenthesis in a string.
		 * 
		 * @param withParen is the String containing the open parenthesis in question.
		 * @param open      is the index of the open parenthesis
		 * 
		 * @return the index of the corresponding close parenthesis
		 * 
		 * @throws ParseException the parse exception
		 */
		private static int findClosingParenthesis(String withParen, int open)
			{
				try
					{
						Stack<Integer> s = new Stack<Integer>();
						for (int i = open + 1; i < withParen.length(); i++)
							{
								char c = withParen.charAt(i);
								switch (c)
									{
									case ')':
										if (s.isEmpty())
											{
												return i;
											} else
											{
												s.pop();
											}
										break;
									case '(':
										s.push(i);
										break;
									}
							}
						throw new ParseException("Missing a close parenthesis.", open);

					} catch (ParseException e)
					{
						e.getStackTrace();
					}
				return -1; // need to allow compile
			}

		/**
		 * Parses the input.
		 * 
		 * @param s is the input string
		 * 
		 * @return a Vector<Object> of all the parsed data as Java Objects
		 * 
		 * @throws IOException Signals that an I/O exception has occurred.
		 */
		public Vector<Object> parse(String s)
			{
				StringTokenizer st = new StringTokenizer(s, "\n");
				final Vector<Object> returning = new Vector<Object>(units.size());
				StringReader sr = new StringReader(st.nextToken());
				for (final Unit u : units)
					{
						final char[] chars = new char[u.length];
						try
							{
								sr.read(chars, 0, u.length);
							} catch (IOException e)
							{
								System.out.println("IOException thrown: Could not read the PDB file. Program terminating.");
								e.printStackTrace();
								System.exit(0);
							}
						StringBuffer sb = new StringBuffer(chars.length);
						for (int i = 0; i < chars.length; i++)
							{
								if (chars[i] != ' ' && chars[i] != 0)
									{
										sb.append(chars[i]);
									}
							}
						String complete = sb.toString();
						if (u.type.equals(R_CHARACTER))
							{
								returning.add(complete);
							} else if (u.type.equals(R_INTEGER))
							{
								returning.add(complete.length() == 0 ? null : Integer.parseInt(complete));
							} else if (u.type.equals(R_REAL_DECIMAL) || u.type.equals(R_REAL_EXPONENT) || u.type.equals(R_REAL_SCIENTIFIC) || u.type.equals(R_REAL_ENGINEERING))
							{
								if (complete.indexOf('E') != -1)
									{
										String end = complete.substring(complete.indexOf("E") + 1);
										if (end.startsWith("+"))
											{
												end = end.substring(1);
											}
										complete = complete.substring(0, complete.indexOf("E"));
										returning.add(complete.length() == 0 ? null
												: Double.parseDouble(complete) / (complete.indexOf('.') == -1 ? Math.pow(10, u.decimalLength) : 1)
														* Math.pow(10, Integer.parseInt(end)));
									} else
									{
										returning.add(
												complete.length() == 0 ? null : Double.parseDouble(complete) / (complete.indexOf('.') == -1 ? Math.pow(10, u.decimalLength) : 1));
									}
							} else if (u.type.equals(R_LOGICAL))
							{
								returning.add(complete.length() == 0 ? null : complete.startsWith("T") || complete.startsWith("t"));
							} else if (u.type.equals(NR_POSITIONING_HORIZONTAL))
							{
								// do nothing
							} else if (u.type.equals(NR_POSITIONING_TAB) || u.type.equals(NR_POSITIONING_TAB_LEFT) || u.type.equals(NR_POSITIONING_TAB_RIGHT))
							{
								// not supported
							} else if (u.type.equals(NR_POSITIONING_VERTICAL))
							{
								sr = new StringReader(st.nextToken());
							}
					}
				return returning;
			}

		/**
		 * Formats the given objects.
		 * 
		 * @param objects are the Java Objects to be formatted
		 * 
		 * @return the formatted string
		 */
		public String format(Vector<Object> objects)
			{
				int minus = 0;
				StringBuffer sb = new StringBuffer();
				int place = -1;
				StringBuffer save = null;
				@SuppressWarnings("unused")
				String sign = "S";
				for (int i = 0; i < objects.size() + minus; i++)
					{
						final Unit u = units.get(i);
						final Object o = objects.get(i - minus);
						if (u.type.equals(R_CHARACTER))
							{
								sb.append(format(o == null ? null : u.length > 0 && ((String) o).length() > u.length ? ((String) o).substring(0, u.length) : (String) o,
										o != null && u.length == 0 ? ((String) o).length() : u.length, true));
							} else if (u.type.equals(R_INTEGER))
							{
								String s = o == null ? null : Integer.toString((Integer) o);
								if (s != null && u.decimalLength > 0)
									{
										boolean neg = s.startsWith("-");
										if (neg)
											{
												s = s.substring(1);
											}
										int numzeros = u.decimalLength - s.length();
										StringBuffer sb2 = new StringBuffer();
										if (neg)
											{
												sb2.append('-');
											}
										for (int j = 0; j < numzeros; j++)
											{
												sb2.append("0");
											}
										sb2.append(s);
										s = sb2.toString();
									}
								sb.append(format(s, u.length, true));
							} else if (u.type.equals(R_REAL_DECIMAL))
							{
								String s = null;
								if (o != null)
									{
										Double d = o instanceof Double ? (Double) o : (Float) o;
										boolean neg = d < 0;
										if (neg)
											{
												d *= -1;
											}
										StringBuffer dfs = new StringBuffer();
										int intLength = Integer.toString(d.intValue()).length();
										for (int j = 0; j < intLength; j++)
											{
												dfs.append('0');
											}
										dfs.append('.');
										for (int j = 0; j < u.decimalLength; j++)
											{
												dfs.append('0');
											}
										s = (neg ? "-" : "") + new DecimalFormat(dfs.toString()).format(d);
									}
								sb.append(format(s, u.length, true));
							} else if (u.type.equals(R_REAL_EXPONENT))
							{
								String s = null;
								if (o != null)
									{
										Double d = o instanceof Double ? (Double) o : (Float) o;
										int exp = 0;
										boolean neg = d < 0;
										if (neg)
											{
												d *= -1;
											}
										while (d > 1)
											{
												d /= 10;
												exp += 1;
											}
										while (d < .1)
											{
												d *= 10;
												exp -= 1;
											}
										boolean expneg = exp < 0;
										if (expneg)
											{
												exp *= -1;
											}
										StringBuffer dfs = new StringBuffer();
										dfs.append("0.");
										for (int j = 0; j < u.decimalLength; j++)
											{
												dfs.append('0');
											}
										s = (neg ? "-" : "") + new DecimalFormat(dfs.toString()).format(d);
										dfs = new StringBuffer();
										for (int j = 0; j < u.exponentLength; j++)
											{
												dfs.append('0');
											}
										s = s + "E" + (expneg ? "-" : "+") + new DecimalFormat(dfs.toString()).format(exp);
									}
								sb.append(format(s, u.length, true));
							} else if (u.type.equals(R_REAL_SCIENTIFIC))
							{
								String s = null;
								if (o != null)
									{
										Double d = o instanceof Double ? (Double) o : (Float) o;
										int exp = 0;
										boolean neg = d < 0;
										if (neg)
											{
												d *= -1;
											}
										while (d > 10)
											{
												d /= 10;
												exp += 1;
											}
										while (d < 1)
											{
												d *= 10;
												exp -= 1;
											}
										boolean expneg = exp < 0;
										if (expneg)
											{
												exp *= -1;
											}
										StringBuffer dfs = new StringBuffer();
										dfs.append("0.");
										for (int j = 0; j < u.decimalLength; j++)
											{
												dfs.append('0');
											}
										s = (neg ? "-" : "") + new DecimalFormat(dfs.toString()).format(d);
										dfs = new StringBuffer();
										for (int j = 0; j < u.exponentLength; j++)
											{
												dfs.append('0');
											}
										s = s + "E" + (expneg ? "-" : "+") + new DecimalFormat(dfs.toString()).format(exp);
									}
								sb.append(format(s, u.length, true));
							} else if (u.type.equals(R_REAL_ENGINEERING))
							{
								String s = null;
								if (o != null)
									{
										Double d = o instanceof Double ? (Double) o : (Float) o;
										int exp = 0;
										boolean neg = d < 0;
										if (neg)
											{
												d *= -1;
											}
										while (d > 10)
											{
												d /= 10;
												exp += 1;
											}
										while (d < 1)
											{
												d *= 10;
												exp -= 1;
											}
										while (exp % 3 != 0)
											{
												d *= 10;
												exp -= 1;
											}
										boolean expneg = exp < 0;
										if (expneg)
											{
												exp *= -1;
											}
										StringBuffer dfs = new StringBuffer();
										dfs.append("0.");
										for (int j = 0; j < u.decimalLength; j++)
											{
												dfs.append('0');
											}
										s = (neg ? "-" : "") + new DecimalFormat(dfs.toString()).format(d);
										dfs = new StringBuffer();
										for (int j = 0; j < u.exponentLength; j++)
											{
												dfs.append('0');
											}
										s = s + "E" + (expneg ? "-" : "+") + new DecimalFormat(dfs.toString()).format(exp);
									}
								sb.append(format(s, u.length, true));
							} else if (u.type.equals(R_LOGICAL))
							{
								String s = o == null ? null : (Boolean) o ? "T" : "F";
								if (s != null)
									{
										StringBuffer sb2 = new StringBuffer();
										for (int j = 0; j < u.length - 1; j++)
											{
												sb2.append(' ');
											}
										sb2.append(s);
										s = sb2.toString();
									}
								sb.append(format(s, u.length, false));
							} else if (u.type.equals(NR_POSITIONING_HORIZONTAL))
							{
								minus++;
								for (int j = 0; j < u.length; j++)
									{
										sb.append(positioningChar);
									}
							} else if (u.type.equals(NR_POSITIONING_TAB) || u.type.equals(NR_POSITIONING_TAB_LEFT) || u.type.equals(NR_POSITIONING_TAB_RIGHT))
							{
								if (save != null)
									{
										while (place - 1 + sb.length() > save.length())
											{
												save.append(' ');
											}
										save.replace(place - 1, place - 1 + sb.length(), sb.toString());
									} else
									{
										save = sb;
									}
								if (u.type.equals(NR_POSITIONING_TAB))
									{
										place = u.length;
									} else if (u.type.equals(NR_POSITIONING_TAB_LEFT))
									{
										place -= u.length - sb.length();
									} else if (u.type.equals(NR_POSITIONING_TAB_RIGHT))
									{
										place += u.length + sb.length();
									}
								sb = new StringBuffer();
								minus++;
							} else if (u.type.equals(NR_POSITIONING_VERTICAL))
							{
								minus++;
								sb.append("\n");
							} else if (u.type.equals(NR_FORMAT_SCANNING_CONTROL))
							{
								minus++;
							} else if (u.type.equals(NR_SIGN_CONTROL_COMPILER) || u.type.equals(NR_SIGN_CONTROL_POSITIVE_ALWAYS) || u.type.equals(NR_SIGN_CONTROL_POSITIVE_NEVER))
							{
								minus++;
								sign = u.type;
								// to vague to implement
							} else if (u.type.equals(NR_BLANK_CONTROL_ZEROS) || u.type.equals(NR_BLANK_CONTROL_REMOVE))
							{
								minus++;
								// not supported
							}
					}
				if (save != null)
					{
						while (place - 1 + sb.length() > save.length())
							{
								save.append(' ');
							}
						save.replace(place - 1, place - 1 + sb.length(), sb.toString());
						sb = save;
						save = null;
						place = -1;
					}
				if (addReturn)
					{
						sb.append("\n");
					}
				return sb.toString();
			}

		/**
		 * Helper method to add spaces and right align content.
		 * 
		 * @param s            is the String to append
		 * @param length       is the desired length
		 * @param rightAligned specifies if the content should be right-aligned
		 * 
		 * @return the formatted string
		 */
		private String format(String s, int length, boolean rightAligned)
			{
				final StringBuffer sb = new StringBuffer();
				if (s == null)
					{
						for (int i = 0; i < length; i++)
							{
								sb.append(' ');
							}
					} else if (length == -1)
					{
						sb.append(s);
					} else if (s.length() > length)
					{
						for (int i = 0; i < length; i++)
							{
								sb.append('*');
							}
					} else
					{
						final int dif = length - s.length();
						if (rightAligned)
							{
								for (int j = 0; j < dif; j++)
									{
										sb.append(" ");
									}
							}
						sb.append(s);
						if (!rightAligned)
							{
								for (int j = 0; j < dif; j++)
									{
										sb.append(' ');
									}
							}
					}
				return sb.toString();
			}

		/**
		 * Gets the positioning char. This is the character to use when skipping spaces during write.
		 * 
		 * @return the positioning char
		 */
		public char getPositioningChar()
			{
				return positioningChar;
			}

		/**
		 * Sets the positioning char. This is the character to use when skipping spaces during write.
		 * 
		 * @param positioningChar the new positioning character
		 */
		public void setPositioningChar(char positioningChar)
			{
				this.positioningChar = positioningChar;
			}

		/**
		 * Specifies whether or not to add a new line at the end of a line during write.
		 * 
		 * @param addReturn the new return line behavior
		 */
		public void setAddReturn(boolean addReturn)
			{
				this.addReturn = addReturn;
			}

		/**
		 * Gets the original format specification string.
		 * 
		 * @return the format specification string
		 */
		public String getFormat()
			{
				return format;
			}

		/**
		 * Checks if returns are added at the end of lines during write.
		 * 
		 * @return true, if is if new lines are added at the end of lines during write
		 */
		public boolean isAddReturn()
			{
				return addReturn;
			}

		/**
		 * Tests many cases and prints the output to standard out.
		 * 
		 * @throws ParseException the parse exception
		 * @throws IOException    Signals that an I/O exception has occurred.
		 */
		public static void test() throws ParseException, IOException
			{
				/*
				 * System.out.println("Fortran Format Tests"); System.out.println();
				 * 
				 * System.out.println( "1. COMMA SANITY CHECKS: (Called during parenthesis expansion. Keep in mind that ambiguous strings may not be interpreted as intended.)");
				 * String check = "I4I4I4"; System.out.println("  a. " + check + " --> " + addCommas(check)); check = "I4I4I4F4.2"; System.out.println("  c. " + check + " --> " +
				 * addCommas(check)); check = "2I4,I4,I4"; System.out.println("  d. " + check + " --> " + addCommas(check)); check = "2I4,5X,I4"; System.out.println("  e. " + check
				 * + " --> " + addCommas(check)); check = "2I4,5X,4I4,2F4.2"; System.out.println("  f. " + check + " --> " + addCommas(check)); check =
				 * "F4.2A5,F4.2A5,2F4.2A5,E4.2E2,4ES5.3E2,4ES5.3E2,E4.2E2"; System.out.println("  g. " + check + " --> " + addCommas(check)); System.out.println();
				 * 
				 * System.out.println("2. PARENTHESIS SANITY CHECKS: (we should see exceptions in some cases)"); check = "I4I4I4"; System.out.print("  a. " + check + " --> "); try
				 * { System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "(I4I4I4)"; System.out.print("  b. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("("
				 * + e.getErrorOffset() + ") " + e.getMessage()); } check = "I4I4I4)"; System.out.print("  c. " + check + " --> "); try {
				 * System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "(I4I4I4"; System.out.print("  d. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" +
				 * e.getErrorOffset() + ") " + e.getMessage()); } check = "2(I4I4I4)"; System.out.print("  e. " + check + " --> "); try {
				 * System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "(I4I4I4)(I4I4I4)"; System.out.print("  f. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) {
				 * System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check = "(I4I4I4(I4I4I4)"; System.out.print("  g. " + check + " --> "); try {
				 * System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "(I4I4I4)I4I4I4)"; System.out.print("  h. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) {
				 * System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check = "(2(I4I4I4)2(A5))"; System.out.print("  i. " + check + " --> "); try {
				 * System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "(2(I4I4I4)A2(A5))"; System.out.print("  j. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) {
				 * System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check = "( 2 ( I 4 I 4 I 4 ) A 2 ( A 5 ) )"; System.out.print("  k. " + check + " --> ");
				 * try { System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "     (I4I4I4)randomtext"; System.out.print("  l. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) {
				 * System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check = "randomtext(I4I4I4)randomtext"; System.out.print("  m. " + check + " --> "); try
				 * { System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "(5(6(A5)))"; System.out.print("  n. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) {
				 * System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check = "(3(5(I2)2(I10)))"; System.out.print("  o. " + check + " --> "); try {
				 * System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "(5(6(A5))"; System.out.print("  p. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("("
				 * + e.getErrorOffset() + ") " + e.getMessage()); } check = "(3(5(I2)2I10)))"; System.out.print("  q. " + check + " --> "); try {
				 * System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check =
				 * "(I4,A2,A4,5X)"; System.out.print("  r. " + check + " --> "); try { System.out.println(removeParenthesis(check)); } catch (ParseException e) {
				 * System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); } check = "(2(I4,3A4,I4),2(A5))"; System.out.print("  s. " + check + " --> "); try {
				 * System.out.println(removeParenthesis(check)); } catch (ParseException e) { System.out.println("(" + e.getErrorOffset() + ") " + e.getMessage()); }
				 * System.out.println();
				 * 
				 * Vector<Object> ints = new Vector<Object>(); ints.add(123); ints.add(-123); ints.add(123456); System.out.println("3. INTEGER (I)");
				 * System.out.println("WRITE: [123,-123,123456]"); String format = "(3I5)"; System.out.println("  a. " + format + " --> '" + write(ints, format) + "'"); format =
				 * "(3I5.2)"; System.out.println("  b. " + format + " --> '" + write(ints, format) + "'"); format = "(3I5.4)"; System.out.println("  c. " + format + " --> '" +
				 * write(ints, format) + "'"); format = "(3I5.5)"; System.out.println("  d. " + format + " --> '" + write(ints, format) + "'"); System.out.println("READ:"); format
				 * = "(4I5)"; String read = "1 3 5 135 135    135"; System.out.println("  a. ('" + read + "') " + format + " --> '" + read(read, format) + "'"); format =
				 * "(I1, I2, I3, I4)"; read = "12 34  56  78  90"; System.out.println("  b. ('" + read + "') " + format + " --> '" + read(read, format) + "'");
				 * System.out.println();
				 * 
				 * Vector<Object> floats = new Vector<Object>(); floats.add(123.345f); floats.add(-123.345f); System.out.println("4. REAL (F)");
				 * System.out.println("WRITE: [123.345,-123.345]"); format = "(2F10.0)"; System.out.println("  a. " + format + " --> '" + write(floats, format) + "'"); format =
				 * "(2F10.1)"; System.out.println("  b. " + format + " --> '" + write(floats, format) + "'"); format = "(2F10.2)"; System.out.println("  c. " + format + " --> '" +
				 * write(floats, format) + "'"); format = "(2F10.3)"; System.out.println("  d. " + format + " --> '" + write(floats, format) + "'"); format = "(2F10.4)";
				 * System.out.println("  e. " + format + " --> '" + write(floats, format) + "'"); format = "(2F10.5)"; System.out.println("  f. " + format + " --> '" +
				 * write(floats, format) + "'"); format = "(2F10.6)"; System.out.println("  g. " + format + " --> '" + write(floats, format) + "'"); format = "(2F10.7)";
				 * System.out.println("  h. " + format + " --> '" + write(floats, format) + "'"); System.out.println("READ:"); format = "(F5.2,F5.2,F5.2)"; read =
				 * "1 2 3 4.5 1 9.4"; System.out.println("  a. ('" + read + "') " + format + " --> '" + read(read, format) + "'"); format = "(F10.4)"; read = "12345E20";
				 * System.out.println("  b. ('" + read + "') " + format + " --> '" + read(read, format) + "'"); format = "(F3.0, F4.1, F6.2, F7.3)"; read = "12 3.4  56E 78.  90";
				 * System.out.println("  c. ('" + read + "') " + format + " --> '" + read(read, format) + "'"); format = "(F6.1, F3.2, F5.0, F6.1)"; read = "12 3.4  56E 78.  90";
				 * System.out.println("  d. ('" + read + "') " + format + " --> '" + read(read, format) + "'"); System.out.println();
				 * 
				 * Vector<Object> doubles = new Vector<Object>(); doubles.add(Math.PI); System.out.println("5. REAL (E)"); System.out.println("WRITE: [\u03C0]"); format =
				 * "(E12.5)"; System.out.println("  a. " + format + " --> '" + write(doubles, format) + "'"); format = "(E12.3E4)"; System.out.println("  b. " + format + " --> '" +
				 * write(doubles, format) + "'"); format = "(E12.7E1)"; System.out.println("  c. " + format + " --> '" + write(doubles, format) + "'"); System.out.println();
				 * 
				 * doubles = new Vector<Object>(); doubles.add(34.5678); System.out.println("6. SCIENTIFIC (ES)"); System.out.println("WRITE: [34.5678]"); format = "(ES12.3E3)";
				 * System.out.println("  a. " + format + " --> '" + write(doubles, format) + "'"); System.out.println();
				 * 
				 * doubles = new Vector<Object>(); doubles.add(1234.567); doubles.add(0.00001234567); System.out.println("7. ENGINEERING (EN)");
				 * System.out.println("WRITE: [1234.567,0.00001234567]"); format = "(2EN12.3E3)"; System.out.println("  a. " + format + " --> '" + write(doubles, format) + "'");
				 * System.out.println();
				 * 
				 * Vector<Object> booleans = new Vector<Object>(); booleans.add(true); booleans.add(false); System.out.println("8. LOGICAL (L)");
				 * System.out.println("WRITE: [true,false]"); format = "(L1,L2)"; System.out.println("  a. " + format + " --> '" + write(booleans, format) + "'"); format =
				 * "(L3,L4)"; System.out.println("  b. " + format + " --> '" + write(booleans, format) + "'"); System.out.println("READ:"); format = "(L3, L8, L10)"; read =
				 * "Fax  Trust   Thursday"; System.out.println("  a. ('" + read + "') " + format + " --> '" + read(read, format) + "'"); System.out.println();
				 * 
				 * Vector<Object> characters = new Vector<Object>(); characters.add("12345"); characters.add("*"); System.out.println("9. CHARACTER (C)");
				 * System.out.println("WRITE: " + characters); format = "(A1,A)"; System.out.println("  a. " + format + " --> '" + write(characters, format) + "'"); format =
				 * "(A2,A)"; System.out.println("  b. " + format + " --> '" + write(characters, format) + "'"); format = "(A3,A)"; System.out.println("  c. " + format + " --> '" +
				 * write(characters, format) + "'"); format = "(A4,A)"; System.out.println("  d. " + format + " --> '" + write(characters, format) + "'"); format = "(A5,A)";
				 * System.out.println("  e. " + format + " --> '" + write(characters, format) + "'"); format = "(A6,A)"; System.out.println("  f. " + format + " --> '" +
				 * write(characters, format) + "'"); format = "(A7,A)"; System.out.println("  g. " + format + " --> '" + write(characters, format) + "'"); format = "(A,A)";
				 * System.out.println("  h. " + format + " --> '" + write(characters, format) + "'"); System.out.println("READ:"); format = "(A4, A5, A7, A)"; read =
				 * "ABCDEFGHIJKLMNOPQRST"; System.out.println("  a. ('" + read + "') " + format + " --> '" + read(read, format) + "'"); format = "(4A5)";
				 * System.out.println("  b. ('" + read + "') " + format + " --> '" + read(read, format) + "'"); System.out.println();
				 * 
				 * Vector<Object> os = new Vector<Object>(); os.add(12); os.add(768); os.add(3.715); System.out.println("10. HORIZONTAL SPACE (X)");
				 * System.out.println("WRITE: [12,768,3.715]"); format = "(1X,I3,3X,I5,2X,F5.2)"; System.out.println("  a. " + format + " --> '" + write(os, format) + "'");
				 * System.out.println();
				 * 
				 * os = new Vector<Object>(); os.add(123); os.add(456); System.out.println("11. TAB (T)"); System.out.println("WRITE: " + os); format = "(T6,I4,T2,I4)";
				 * System.out.println("  a. " + format + " --> '" + write(os, format) + "'"); System.out.println();
				 * 
				 * os = new Vector<Object>(); os.add(123); os.add(456); os.add(789); System.out.println("12. TABS (T,TL,TR)"); System.out.println("WRITE: " + os); format =
				 * "(T10,I3,TL9,I3,TR5,I3)"; System.out.println("  a. " + format + " --> '" + write(os, format) + "'"); System.out.println();
				 * 
				 * os = new Vector<Object>(); os.add(123); os.add(456);
				 */

				// os.add("+-*/");

				/*
				 * System.out.println("13. VERTICAL POSITIONING (/)"); System.out.println("WRITE: " + os); format = "(I5//I6/1X,A)"; System.out.println("  a. " + format + " --> '"
				 * + write(os, format) + "'"); System.out.println("READ:"); format = "(I5/I5/I5)"; read = "  123  456\n  789  012\n  345  678"; System.out.println("  a. ('" + read
				 * + "') " + format + " --> '" + read(read, format) + "'"); System.out.println();
				 */
			}

	}
