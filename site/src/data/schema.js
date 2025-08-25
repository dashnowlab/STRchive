import Ajv from "ajv/dist/2020";
import { customizeValidator } from "@rjsf/validator-ajv8";

export const validator = customizeValidator({
  ajvOptionsOverrides: { strict: false },
  AjvClass: Ajv,
});

/** add custom keyword for validating one value less than other value */
validator.ajv.addKeyword({
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
          params: { deps: null },
        },
      ];
    else validate.errors = null;
    return valid;
  },
  errors: true,
});

/** add custom keyword for validating one value greater than other value */
validator.ajv.addKeyword({
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
          params: { deps: null },
        },
      ];
    else validate.errors = null;
    return valid;
  },
  errors: true,
});
