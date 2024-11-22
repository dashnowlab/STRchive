/** get full name of inheritance from abbreviation */
export const getInheritance = (abbreviation) => {
  if (abbreviation == "AD") return "Autosomal dominant";
  if (abbreviation == "AR") return "Autosomal recessive";
  if (abbreviation == "XLD" || abbreviation == "XD") return "X-linked dominant";
  if (abbreviation == "XLR" || abbreviation == "XR")
    return "X-linked recessive";
  return "";
};
