import Ajv, { ValidationError } from "ajv/dist/2020";
import schema from "~/STRchive-loci.schema.json";

/** validator library */
export const ajv = new Ajv({ strict: false, allErrors: true });

/** add custom keyword for validating one value less than other value */
ajv.addKeyword({
  keyword: "less_than",
  validate: function validate(otherKey, value, parentSchema, { parentData }) {
    /** value of referenced property */
    const otherValue = parentData[otherKey];
    /** if null is allowable and if either value null, don't error */
    if (
      [parentSchema.type].flat()?.includes("null") &&
      (value === null || otherValue === null)
    )
      return true;
    /** test */
    const valid = value <= otherValue;
    /** custom error message */
    if (!valid)
      validate.errors = [
        {
          keyword: "less_than",
          message: `Value must be less than ${otherKey}`,
        },
      ];
    return valid;
  },
  errors: true,
});

/** add custom keyword for validating one value greater than other value */
ajv.addKeyword({
  keyword: "greater_than",
  validate: function validate(otherKey, value, parentSchema, { parentData }) {
    /** value of referenced property */
    const otherValue = parentData[otherKey];
    /** if null is allowable and if either value null, don't error */
    if (
      [parentSchema.type].flat()?.includes("null") &&
      (value === null || otherValue === null)
    )
      return true;
    /** test */
    const valid = value >= otherValue;
    /** custom error message */
    if (!valid)
      validate.errors = [
        {
          keyword: "greater_than",
          message: `Value must be greater than ${otherKey}`,
        },
      ];
    return valid;
  },
  errors: true,
});

/** compile schema */
export const validator = ajv.compile(schema);

/** validator function */
export const validate = (data) => {
  validator(data);
  return validator.errors;
};
