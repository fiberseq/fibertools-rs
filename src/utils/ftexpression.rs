use crate::fiber;
use crate::utils::input_bam;

#[derive(Debug, PartialEq, Clone)]
enum Comparison {
    LessThan,
    LessThanOrEqual,
    GreaterThan,
    GreaterThanOrEqual,
    Equal,
}

#[derive(Debug, Clone)]
enum FunctionArgument {
    LenArgument(String),
    QualArgument(String),
}

#[derive(Debug, Clone)]
enum Expression {
    FunctionCall {
        name: String,
        argument: FunctionArgument,
    },
    Comparison {
        left: Box<Expression>,
        operator: Comparison,
        right: f64,
    },
    And(Box<Expression>, Box<Expression>),
}

impl Expression {
    fn evaluate(&self) -> bool {
        match self {
            Expression::FunctionCall { name: _, argument } => {
                let arg_len = match argument {
                    FunctionArgument::LenArgument(arg) => arg.len(),
                    FunctionArgument::QualArgument(arg) => arg.len(),
                };
                println!("FunctionCall result: {}", arg_len);
                true
            }
            Expression::Comparison {
                left,
                operator,
                right,
            } => {
                let left_value = match &**left {
                    Expression::FunctionCall { name, argument } => match name.as_str() {
                        "len" | "qual" => match argument {
                            FunctionArgument::LenArgument(arg)
                            | FunctionArgument::QualArgument(arg) => arg.len() as f64,
                        },
                        _ => panic!("Unknown function name"),
                    },
                    _ => panic!("Expected function call on left side of comparison"),
                };

                match operator {
                    Comparison::LessThan => left_value < *right,
                    Comparison::LessThanOrEqual => left_value <= *right,
                    Comparison::GreaterThan => left_value > *right,
                    Comparison::GreaterThanOrEqual => left_value >= *right,
                    Comparison::Equal => left_value == *right,
                }
            }
            Expression::And(expr1, expr2) => expr1.evaluate() && expr2.evaluate(),
        }
    }
}

struct Parser {
    input: String,
    pos: usize,
}

impl Parser {
    fn new(input: &str) -> Self {
        Parser {
            input: input.chars().filter(|c| !c.is_whitespace()).collect(),
            pos: 0,
        }
    }

    fn parse(&mut self) -> Result<Expression, String> {
        let function_call = self.parse_function_call()?;
        if self.consume_if("=:") {
            let (min, max) = self.parse_range()?;
            return Ok(Expression::And(
                Box::new(Expression::Comparison {
                    left: Box::new(function_call.clone()),
                    operator: Comparison::GreaterThanOrEqual,
                    right: min,
                }),
                Box::new(Expression::Comparison {
                    left: Box::new(function_call),
                    operator: Comparison::LessThan,
                    right: max,
                }),
            ));
        }
        let operator = self.parse_comparison_operator()?;
        let number = self.parse_number()?;
        Ok(Expression::Comparison {
            left: Box::new(function_call),
            operator,
            right: number,
        })
    }

    fn parse_function_call(&mut self) -> Result<Expression, String> {
        let name = self.parse_identifier()?;
        self.expect_char('(')?;
        let argument = match name.as_str() {
            "len" | "qual" => {
                let arg = self.parse_identifier()?;
                match name.as_str() {
                    "len" if matches!(arg.as_str(), "msp" | "fire" | "nuc" | "lnk") => {
                        FunctionArgument::LenArgument(arg)
                    }
                    "qual" if matches!(arg.as_str(), "m6a" | "5mC") => {
                        FunctionArgument::QualArgument(arg)
                    }
                    _ => return Err(format!("Invalid argument for {} function: {}", name, arg)),
                }
            }
            _ => return Err(format!("Unknown function name: {}", name)),
        };
        self.expect_char(')')?;
        Ok(Expression::FunctionCall { name, argument })
    }

    fn parse_identifier(&mut self) -> Result<String, String> {
        let start = self.pos;
        while self.peek_char().map_or(false, |c| c.is_alphanumeric()) {
            self.pos += 1;
        }
        if start == self.pos {
            Err(format!("Expected identifier at position {}", self.pos))
        } else {
            Ok(self.input[start..self.pos].to_string())
        }
    }

    fn parse_comparison_operator(&mut self) -> Result<Comparison, String> {
        let operator = match self.next_char() {
            Some('<') => match self.peek_char() {
                Some('=') => {
                    self.pos += 1;
                    Comparison::LessThanOrEqual
                }
                _ => Comparison::LessThan,
            },
            Some('>') => match self.peek_char() {
                Some('=') => {
                    self.pos += 1;
                    Comparison::GreaterThanOrEqual
                }
                _ => Comparison::GreaterThan,
            },
            Some('=') => Comparison::Equal,
            _ => {
                return Err(format!(
                    "Expected comparison operator at position {}",
                    self.pos
                ))
            }
        };
        Ok(operator)
    }

    fn parse_number(&mut self) -> Result<f64, String> {
        let start = self.pos;
        while self
            .peek_char()
            .map_or(false, |c| c.is_digit(10) || c == '.')
        {
            self.pos += 1;
        }
        self.input[start..self.pos]
            .parse::<f64>()
            .map_err(|e| format!("Failed to parse number: {}", e))
    }

    fn parse_range(&mut self) -> Result<(f64, f64), String> {
        let min = self.parse_number()?;
        self.expect_char(':')?;
        let max = self.parse_number()?;
        Ok((min, max))
    }

    fn expect_char(&mut self, expected: char) -> Result<(), String> {
        match self.next_char() {
            Some(c) if c == expected => Ok(()),
            _ => Err(format!("Expected '{}' at position {}", expected, self.pos)),
        }
    }

    fn consume_if(&mut self, s: &str) -> bool {
        if self.input[self.pos..].starts_with(s) {
            self.pos += s.len();
            true
        } else {
            false
        }
    }

    fn peek_char(&self) -> Option<char> {
        self.input[self.pos..].chars().next()
    }

    fn next_char(&mut self) -> Option<char> {
        let ch = self.peek_char();
        self.pos += ch.map_or(0, |c| c.len_utf8());
        ch
    }
}

pub fn filter_fiber(fiber: &mut fiber::FiberseqData, input_bam_filter: &input_bam::FiberFilters) {
    let expression = match input_bam_filter.filter_expression {
        Some(ref expression) => expression,
        None => return,
    };
    let mut parser = Parser::new(&expression);
    // ideally I think parser should be filtering one of "fiber"s ranges, msp, nuc, cpg, m6a.
    // parser my need a "name" for this
    todo!()
}

/// @SHANE you can run these tests with
/// cargo test --lib -- utils::ftexpression::tests --show-output
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parser() {
        let expressions = vec![
            "len(msp)>1000",
            "len(fire)>=500",
            "qual(m6a)=0.7",
            "qual(5mC)<1.0",
            "len(msp)=30:40",
        ];

        for expression in expressions {
            let mut parser = Parser::new(expression);
            match parser.parse() {
                Ok(expr) => {
                    println!("Expression: {}", expression);
                    println!("Evaluates to: {}", expr.evaluate());
                }
                Err(e) => println!("Failed to parse expression: {}. Error: {}", expression, e),
            }
        }
    }
}
